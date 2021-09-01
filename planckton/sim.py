"""Tools for running an OPV simulation with PlanckTon."""
import os

import hoomd.data
import hoomd.dump
import hoomd.md
import unyt as u
from mbuild.formats.hoomd_simulation import create_hoomd_simulation

from planckton.utils.solvate import set_coeffs


class Simulation:
    """
    Convenience class for initializing and running a HOOMD simulation.

    Parameters
    ----------
    typed_system : ParmEd structure
        Typed structure used to initialize the simulation
    kT : float
        Dimensionless temperature at which to run the simulation
    e_factor : float, default 1.0
        Scaling parameter for particle interaction strengths, used to simulate
        solvent
    tau : float, default 5.0
        Thermostat coupling period
    r_cut : float, default 2.5
        Cutoff radius for potentials
    gsd_write : int, default 1e6
        Period to write simulation snapshots to gsd file
    log_write : int, default 1e5
        Period to write simulation data to the log file
    shrink_steps : int, default 1e6
        Number of timesteps over which to shrink the box
    shrink_kT_reduced : float, default 10
        Dimensionless temperature to run the shrink step
    n_steps : int, default 1e3
        Number of steps to run the simulation
    dt : float, default 0.0001
        Size of simulation timestep in simulation time units
    mode : str, default "gpu"
        Mode flag passed to hoomd.context.initialize. Options are "cpu" and
        "gpu".
    target_length : unyt.unyt_quantity, default None
        Target final box length for the shrink step. If None is provided, no
        shrink step will be performed.
    restart : str, default None
        Path to gsd file from which to restart the simulation

    Attributes
    ----------
    system : ParmEd structure
        Structure used to initialize the simulation
    kT : float
        Dimensionless temperature at the simulation is run
    e_factor : float
        Scaling parameter for particle interaction strengths
    tau : float
        Thermostat coupling period
    r_cut : float
        Cutoff radius for potentials
    gsd_write : int
        Period to write simulation snapshots to gsd file
    log_write : int
        Period to write simulation data to the log file
    shrink_steps : int
        Number of timesteps over which to shrink the box
    shrink_kT_reduced : float
        Dimensionless temperature to run the shrink step
    n_steps : int
        Number of steps to run the simulation
    dt : float
        Size of simulation timestep in simulation time units
    mode : str
        Mode flag passed to hoomd.context.initialize.
    target_length : unyt.unyt_quantity
        Target final box length for the shrink step.
    """

    def __init__(
        self,
        typed_system,
        kT,
        e_factor=1.0,
        tau=5.0,
        r_cut=2.5,
        gsd_write=1e5,
        log_write=1e3,
        shrink_steps=1e3,
        shrink_kT_reduced=10,
        n_steps=1e7,
        dt=0.0001,
        mode="gpu",
        target_length=None,
        restart=None,
    ):
        self.system = typed_system
        self.kT = kT
        self.e_factor = e_factor
        self.tau = tau
        self.r_cut = r_cut
        self.gsd_write = gsd_write
        self.log_write = log_write
        self.shrink_steps = shrink_steps
        self.shrink_kT_reduced = shrink_kT_reduced
        self.n_steps = n_steps
        self.dt = dt
        self.mode = mode
        self.target_length = target_length
        self.restart = restart

    def run(self):
        """Run the simulation."""
        hoomd_args = f"--single-mpi --mode={self.mode}"
        sim = hoomd.context.initialize(hoomd_args)

        with sim:
            hoomd.util.quiet_status()
            # mbuild units are nm, amu
            hoomd_objects, ref_values = create_hoomd_simulation(
                self.system,
                auto_scale=True,
                restart=self.restart,
                r_cut=self.r_cut,
            )
            self.ref_values = ref_values
            snap = hoomd_objects[0]
            hoomd.util.unquiet_status()

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
                    group=both, tau=self.tau, kT=self.shrink_kT_reduced
                )
            except NameError:
                # both does not exist
                integrator = hoomd.md.integrate.nvt(
                    group=all_particles, tau=self.tau, kT=self.shrink_kT_reduced
                )

            hoomd.dump.gsd(
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
            if self.restart is None:
                integrator.randomize_velocities(seed=42)

            if self.target_length is not None:
                # Run the shrink step
                final_length = self.target_length.to("Angstrom").value
                final_box = (self.shrink_steps, final_length)
                size_variant = hoomd.variant.linear_interp(
                    [(0, snap.box.Lx), final_box], zero=0
                )
                box_resize = hoomd.update.box_resize(L=size_variant)
                hoomd.run_upto(self.shrink_steps)
                box_resize.disable()
                self.n_steps += self.shrink_steps

            # After shrinking, reset velocities and change temp
            integrator.set_params(kT=self.kT)
            integrator.randomize_velocities(seed=42)
            integrator_mode.set_params(dt=self.dt)

            try:
                hoomd.run_upto(self.n_steps + 1, limit_multiple=self.gsd_write)
                print("Simulation completed")
                done = True
            except hoomd.WalltimeLimitReached:
                print("Walltime limit reached")
                done = False
            finally:
                gsd_restart.write_restart()
                print("Restart file written")
            return done
