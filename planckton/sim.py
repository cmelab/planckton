"""Tools for running an OPV simulation with PlanckTon."""
import os

import hoomd
import unyt as u
from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield


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
    shrink_period: int, default 1
        Number of steps to run in between box updates.
    target_length : unyt.unyt_quantity, default None
        Target final box length for the shrink step. If None is provided, no
        shrink step will be performed.
    restart : str, default None
        Path to gsd file from which to restart the simulation.
    nlist : str, default "Cell"
        Type of neighborlist to use. Options are "Cell", "Tree", and Stencil".
    seed: int, default 5
        Random seed for hoomd simulation

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
    e_factor : float, default 1.0
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
    shrink_period: int, default 1
        Number of steps to run in between box updates.
    target_length : unyt.unyt_quantity
        Target final box length for the shrink step.
    nlist : hoomd.md.nlist
        Type of neighborlist used, see
        https://hoomd-blue.readthedocs.io/en/stable/module-md-nlist.html
        for more information.
    seed: int, default 5
        Random seed for hoomd simulation
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
        shrink_period=1,
        target_length=None,
        restart=None,
        nlist="Cell",
        seed=5,
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
        self.shrink_period = shrink_period
        self.n_steps = n_steps
        self.dt = dt
        if target_length:
            self.target_length = target_length.to("Angstrom").value
        else:
            self.target_length = None
        self.restart = restart
        self.nlist = getattr(hoomd.md.nlist, nlist)
        self.seed = seed
        self.log_quantities = [
            "kinetic_temperature",
            "pressure",
            "volume",
            "potential_energy",
            "kinetic_energy",
        ]

    def run(self):
        """Run the simulation."""
        device = hoomd.device.auto_select()
        sim = hoomd.Simulation(device=device, seed=self.seed)

        # mbuild units are nm, amu
        snap, hoomd_forcefield, ref_values = create_hoomd_forcefield(
            self.system, auto_scale=True, r_cut=self.r_cut
        )

        lj = [
            force
            for force in hoomd_forcefield
            if isinstance(force, hoomd.md.pair.LJ)
        ][0]
        if self.e_factor != 1:
            print(f"Scaling LJ epsilon values for all pairs by {self.e_factor}")
            for pair in lj.params:
                lj.params[pair]["epsilon"] *= self.e_factor

        if not isinstance(self.nlist, hoomd.md.nlist.Cell):
            exclusions = lj.nlist.exclusions
            lj.nlist = self.nlist(buffer=0.4)
            lj.nlist.exclusions = exclusions

        if self.restart:
            sim.create_state_from_gsd(self.restart)
        else:
            sim.create_state_from_snapshot(snap)

        if self.target_length:
            self.target_length /= ref_values.distance

        integrator = hoomd.md.Integrator(dt=self.dt)
        all_particles = hoomd.filter.All()
        integrator_method = hoomd.md.methods.NVT(
            filter=all_particles, kT=self.shrink_kT, tau=self.shrink_tau
        )
        integrator.forces = hoomd_forcefield
        integrator.methods = [integrator_method]
        sim.operations.add(integrator)

        gsd_writer, table_file = self._hoomd_writers(
            group=all_particles, sim=sim, forcefields=hoomd_forcefield
        )
        sim.operations.writers.append(gsd_writer)
        sim.operations.writers.append(table_file)

        if self.target_length is not None:
            # Run the shrink step
            box_resize_trigger = hoomd.trigger.Periodic(self.shrink_period)
            ramp = hoomd.variant.Ramp(
                A=0, B=1, t_start=0, t_ramp=int(self.shrink_steps)
            )
            initial_box = sim.state.box
            final_box = hoomd.Box(
                Lx=self.target_length,
                Ly=self.target_length,
                Lz=self.target_length,
            )
            box_resize = hoomd.update.BoxResize(
                box1=initial_box,
                box2=final_box,
                variant=ramp,
                trigger=box_resize_trigger,
            )
            sim.operations.updaters.append(box_resize)
            sim.state.thermalize_particle_momenta(
                filter=all_particles, kT=self.shrink_kT
            )
            sim.run(self.shrink_steps, write_at_start=True)

            # Begin temp ramp
            for kT, tau, n_steps in zip(self.kT, self.tau, self.n_steps):
                sim.operations.integrator.methods[0].kT = kT
                sim.operations.integrator.methods[0].tau = tau
                sim.state.thermalize_particle_momenta(
                    filter=all_particles, kT=kT
                )

                sim.run(n_steps)
                if sim.timestep >= sum(self.n_steps) + self.shrink_steps:
                    print("Simulation completed")
                    done = True
                else:
                    print("Simulation not completed.")
                    done = False
                hoomd.write.GSD.write(
                    state=sim.state, mode="wb", filename="restart.gsd"
                )
            return done

    def _hoomd_writers(self, group, forcefields, sim):
        # GSD and Logging:
        if self.restart:
            writemode = "a"
        else:
            writemode = "w"
        gsd_writer = hoomd.write.GSD(
            filename="trajectory.gsd",
            trigger=hoomd.trigger.Periodic(period=int(self.gsd_write), phase=0),
            mode=f"{writemode}b",
            dynamic=["momentum"],
        )
        logger = hoomd.logging.Logger(categories=["scalar", "string"])
        logger.add(sim, quantities=["timestep", "tps"])
        thermo_props = hoomd.md.compute.ThermodynamicQuantities(filter=group)
        sim.operations.computes.append(thermo_props)
        logger.add(thermo_props, quantities=self.log_quantities)
        for f in forcefields:
            logger.add(f, quantities=["energy"])

        table_file = hoomd.write.Table(
            output=open("trajectory.txt", mode=f"{writemode}", newline="\n"),
            trigger=hoomd.trigger.Periodic(period=int(self.log_write), phase=0),
            logger=logger,
            max_header_len=None,
        )
        return gsd_writer, table_file
