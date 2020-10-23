import logging
import os

import hoomd.data
import hoomd.deprecated
import hoomd.dump
import hoomd.md
from cme_utils.manip.convert_rigid import init_wrapper
from cme_utils.manip.ff_from_foyer import set_coeffs


class Simulation:
    """
    Convenience class for initializing and running a HOOMD simulation.

    Parameters
    ----------
    input_xml : str
        Hoomdxml file to use to initialize the simulation
    kT : float
        Dimensionless temperature at which to run the simulation
    e_factor : float
        Scaling parameter for particle interaction strengths, used to simulate
        solvent (default 1.0)
    tau : float
        Thermostat coupling strength (default 5.0)
    gsd_write : int
        Period to write simulation snapshots to gsd file (default 1e6)
    log_write : int
        Period to write simulation data to the log file (default 1e5)
    shrink_time : int
        Number of timesteps over which to shrink the box (default 1e6)
    shrink_kT_reduced : float
        Dimensionless temperature to run the shrink step (default 10)
    n_steps : int
        Number of steps to run the simulation (default 1e3)
    dt : float
        Size of simulation timestep in simulation time units (default 0.0001)
    mode : str
        Mode flag passed to hoomd.context.initialize. Options are "cpu" and
        "gpu". (default "gpu")
    target_length : float
        Target final box length for the shrink step. If None is provided, no
        shrink step will be performed. (default None)

    Attributes
    ----------
    input_xml : str
        Hoomdxml file used to initialize the simulation
    kT : float
        Dimensionless temperature at the simulation is run
    e_factor : float
        Scaling parameter for particle interaction strengths
    tau : float
        Thermostat coupling strength
    gsd_write : int
        Period to write simulation snapshots to gsd file
    log_write : int
        Period to write simulation data to the log file
    shrink_time : int
        Number of timesteps over which to shrink the box
    shrink_kT_reduced : float
        Dimensionless temperature to run the shrink step
    n_steps : int
        Number of steps to run the simulation
    dt : float
        Size of simulation timestep in simulation time units
    mode : str
        Mode flag passed to hoomd.context.initialize.
    target_length : float
        Target final box length for the shrink step.
    """
    def __init__(
        self,
        input_xml,
        kT,
        e_factor=1.0,
        tau=5.0,
        gsd_write=1e6,
        log_write=1e5,
        shrink_time=1e6,
        shrink_kT_reduced=10,
        n_steps=1e3,
        dt=0.0001,
        mode="gpu",
        target_length=None,
    ):
        self.input_xml = input_xml
        self.e_factor = e_factor
        self.tau = tau
        self.kT = kT
        self.gsd_write = gsd_write
        self.log_write = log_write
        self.shrink_time = shrink_time
        self.shrink_kT_reduced = shrink_kT_reduced
        self.n_steps = n_steps
        self.dt = dt
        self.mode = mode
        self.target_length = target_length

    def run(self):
        if hoomd.context.exec_conf is None:
            hoomd_args = f"--single-mpi --mode={self.mode}"
            hoomd.context.initialize(hoomd_args)
        with hoomd.context.SimulationContext():
            # TODO Robust restart logic when reading in rigid bodies
            if os.path.isfile("restart.gsd"):
                system = hoomd.init.read_gsd(
                        filename=None, restart="restart.gsd"
                        )
            else:
                system = init_wrapper(self.input_xml)
            nl = hoomd.md.nlist.cell()
            logging.info("Setting coeffs")
            hoomd.util.quiet_status()
            system = set_coeffs(self.input_xml, system, nl, self.e_factor)
            hoomd.util.unquiet_status()
            integrator_mode = hoomd.md.integrate.mode_standard(dt=self.dt)
            rigid = hoomd.group.rigid_center()
            nonrigid = hoomd.group.nonrigid()
            both_group = hoomd.group.union("both", rigid, nonrigid)
            all_particles = hoomd.group.all()
            integrator = hoomd.md.integrate.nvt(
                group=both_group, tau=self.tau, kT=self.shrink_kT_reduced
            )
            hoomd.dump.gsd(
                filename="trajectory.gsd",
                period=self.gsd_write,
                group=all_particles,
                overwrite=False,
                phase=0,
            )
            gsd_restart = hoomd.dump.gsd(
                "restart.gsd",
                period=self.gsd_write,
                group=all_particles,
                truncate=True,
                phase=0,
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
            integrator.randomize_velocities(seed=42)

            if self.target_length == None:
                self.target_length = system.box.Lx
            size_variant = hoomd.variant.linear_interp(
                [(0, system.box.Lx), (self.shrink_time, self.target_length)], zero=0
            )
            box_resize = hoomd.update.box_resize(L=size_variant)
            hoomd.run_upto(self.shrink_time)
            box_resize.disable()

            # After shrinking, reset velocities and change temp
            integrator.set_params(kT=self.kT)
            integrator.randomize_velocities(seed=42)
            integrator_mode.set_params(dt=self.dt)

            try:
                hoomd.run_upto(self.n_steps + 1, limit_multiple=self.gsd_write)
            except hoomd.WalltimeLimitReached:
                pass
            finally:
                gsd_restart.write_restart()
                hoomd.deprecated.dump.xml(
                    group=hoomd.group.all(), filename="final.xml", all=True
                )
