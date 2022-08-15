"""Tools for running an OPV simulation with PlanckTon."""
import os

import hoomd.data
# import hoomd.dump
import hoomd.md
import unyt as u
from mbuild.formats.hoomd_forcefield import create_hoomd_ff

from planckton.utils.solvate import set_coeffs


class Simulation:
    """Convenience class for initializing and running a HOOMD simulation.

    Parameters
    ----------
    typed_system : ParmEd structure
        Typed structure used to initialize the simulation.
    kT : list of float
        Dimensionless temperature(s) at which to run the simulation
    tau : list of float
        Thermostat coupling period(s) (in simulation time units)
    n_steps : list of int
        Number of timesteps to run each simulation block
    dt : float, default 0.001
        Size of simulation timestep (in simulation time units)
    e_factor : float, default 1.0
        Scaling parameter for particle interaction strengths, used to simulate
        solvent
    r_cut : float, default 2.5
        Cutoff radius for potentials (in simulation distance units)
    gsd_write : int, default 1e6
        Period to write simulation snapshots to gsd file
    log_write : int, default 1e5
        Period to write simulation data to the log file
    shrink_steps : int, default 1e6
        Number of timesteps over which to shrink the box
    shrink_kT : float, default 10
        Dimensionless temperature to run the shrink step
        "gpu".
    target_length : unyt.unyt_quantity, default None
        Target final box length for the shrink step. If None is provided, no
        shrink step will be performed.
    restart : str, default None
        Path to gsd file from which to restart the simulation.
    nlist : str, default "cell"
        Type of neighborlist to use. Options are "cell", "tree", and "stencil".
        See https://hoomd-blue.readthedocs.io/en/stable/nlist.html and
        https://hoomd-blue.readthedocs.io/en/stable/module-md-nlist.html

    Attributes
    ----------
    ref_values : namedtuple
        Distance, energy, and mass values used for scaling in angstroms,
        kcal/mol, and daltons.
    system : ParmEd structure
        Structure used to initialize the simulation
    kT : list of float
        Dimensionless temperature(s) at the simulation is run
    tau : list of float
        Thermostat coupling period(s)
    n_steps : list of int
        Number of timesteps to run each simulation block
    dt : float
        Size of simulation timestep in simulation time units
    e_factor : float
        Scaling parameter for particle interaction strengths
    r_cut : float
        Cutoff radius for potentials
    gsd_write : int
        Period to write simulation snapshots to gsd file
    log_write : int
        Period to write simulation data to the log file
    shrink_steps : int
        Number of timesteps over which to shrink the box
    shrink_kT : float
        Dimensionless temperature to run the shrink step
    shrink_tau : float
        Thermostat coupling period during shrink step
    target_length : unyt.unyt_quantity
        Target final box length for the shrink step.
    nlist : hoomd.md.nlist
        Type of neighborlist used, see
        https://hoomd-blue.readthedocs.io/en/stable/module-md-nlist.html
        for more information.
    """

    def __init__(
        self,
        typed_system,
        kT,
        tau,
        n_steps,
        dt=0.001,
        e_factor=1.0,
        r_cut=2.5,
        gsd_write=1e5,
        log_write=1e3,
        shrink_steps=1e3,
        shrink_kT=10,
        shrink_tau=1.0,
        target_length=None,
        restart=None,
        nlist="Cell",
    ):
        assert len(kT) == len(tau) == len(n_steps), (
            f"Must have the same number of values for kT (found {len(kT)}), "
            f"tau (found {len(tau)}), and n_steps (found {len(n_steps)})."
        )

        # Combine n_steps, so each value reflects the total steps at that point
        for i, n in enumerate(n_steps[:-1]):
            for j in range(i + 1, len(n_steps)):
                n_steps[j] += n

        self.system = typed_system
        self.kT = kT
        self.e_factor = e_factor
        self.tau = tau
        self.r_cut = r_cut
        self.gsd_write = gsd_write
        self.log_write = log_write
        self.shrink_steps = shrink_steps
        self.shrink_kT = shrink_kT
        self.shrink_tau = shrink_tau
        self.n_steps = n_steps
        self.dt = dt
        self.target_length = target_length
        self.restart = restart
        self.nlist = getattr(hoomd.md.nlist, nlist)

    def run(self):
        """Run the simulation."""
        device = hoomd.device.auto_select()
        sim = hoomd.Simulation(device=device)

        # mbuild units are nm, amu
        snap, hoomd_objects, ref_values = create_hoomd_ff(
            self.system,
            auto_scale=True,
            restart=self.restart,
            nlist=self.nlist,
            r_cut=self.r_cut,
        )
        snap = hoomd_objects[0]

        if self.target_length is not None:
            self.target_length /= ref_values.distance

        if self.e_factor != 1:
            print("Scaling LJ coeffs by e_factor")
            hoomd.util.quiet_status()
            # catch all instances of LJ pair
            ljtypes = [
                i
                for i in sim.forces
                if isinstance(i, hoomd.md.pair.lj)
                or isinstance(i, hoomd.md.special_pair.lj)
            ]

            for lj in ljtypes:
                pair_list = lj.get_metadata()["pair_coeff"].get_metadata()
                for pair_dict in pair_list:
                    # Scale the epsilon values by e_factor
                    try:
                        a, b, new_dict = set_coeffs(
                            pair_dict, self.e_factor
                        )
                        lj.pair_coeff.set(a, b, **new_dict)
                    except ValueError:
                        # if the pair has not been defined,
                        # it will not have a dictionary object
                        # instead it will be a string (e.g. "ca-ca")
                        # and will fail when trying to make the new_dict
                        pass
            hoomd.util.unquiet_status()

            integrator_mode = hoomd.md.integrate.mode_standard(dt=self.dt)
            all_particles = hoomd.group.all()
            try:
                integrator = hoomd.md.integrate.nvt(
                    group=both, tau=self.shrink_tau, kT=self.shrink_kT
                )
            except NameError:
                # both does not exist
                integrator = hoomd.md.integrate.nvt(
                    group=all_particles, tau=self.shrink_tau, kT=self.shrink_kT
                )

            hoomd.writer.gsd(
                filename="trajectory.gsd",
                period=self.gsd_write,
                group=all_particles,
                overwrite=False,
                phase=0,
                dynamic=["momentum"],
            )
            gsd_restart = hoomd.dump.gsd(
                "restart.gsd",
                period=self.gsd_write,
                group=all_particles,
                truncate=True,
                phase=0,
                dynamic=["momentum"],
            )
            log_quantities = [
                "temperature",
                "pressure",
                "volume",
                "potential_energy",
                "kinetic_energy",
                "pair_lj_energy",
                "bond_harmonic_energy",
                "angle_harmonic_energy",
            ]
            hoomd.analyze.log(
                "trajectory.log",
                quantities=log_quantities,
                period=self.log_write,
                header_prefix="#",
                overwrite=False,
                phase=0,
            )

            if self.target_length is not None:
                # Run the shrink step
                final_length = self.target_length.to("Angstrom").value
                final_box = (self.shrink_steps, final_length)
                size_variant = hoomd.variant.linear_interp(
                    [(0, snap.box.Lx), final_box], zero=0
                )
                box_resize = hoomd.update.box_resize(L=size_variant)
                integrator.randomize_velocities(seed=42)
                hoomd.run_upto(self.shrink_steps)
                box_resize.disable()
                self.n_steps = [i + self.shrink_steps for i in self.n_steps]

            # Begin temp ramp
            for kT, tau, n_steps in zip(self.kT, self.tau, self.n_steps):
                integrator.set_params(kT=kT, tau=tau)
                # Reset velocities
                integrator.randomize_velocities(seed=42)

                try:
                    hoomd.run_upto(n_steps + 1, limit_multiple=self.gsd_write)
                    if sim.system.getCurrentTimeStep() >= self.n_steps[-1]:
                        print("Simulation completed")
                        done = True
                except hoomd.WalltimeLimitReached:
                    print("Walltime limit reached")
                    done = False
                finally:
                    gsd_restart.write_restart()
                    print("Restart file written")
        return done
