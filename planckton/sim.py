import logging
import os

import hoomd.data
import hoomd.deprecated
import hoomd.dump
import hoomd.md
from cme_utils.manip.convert_rigid import init_wrapper
from cme_utils.manip.ff_from_foyer import set_coeffs


class Simulation:
    def __init__(
        self,
        kT,
        e_factor=1.0,
        tau=5.0,
        gsd_write=1e6,
        log_write=1e5,
        shrink_time=1e6,
        shrink_factor=5,
        shrink_kT_reduced=10,
        n_steps=1e3,
        dt=0.0001,
        mode="gpu",
        target_length=None,
    ):
        self.e_factor = e_factor
        self.tau = tau
        self.kT = kT
        self.gsd_write = gsd_write
        self.log_write = log_write
        self.shrink_time = shrink_time
        self.shrink_factor = shrink_factor
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
                system = hoomd.init.read_gsd(filename=None, restart="restart.gsd")
            else:
                system = init_wrapper(self.input_xml)
            nl = hoomd.md.nlist.cell()
            logging.info("Setting coefs")
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
                "dihedral_table_energy",
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
                [(0, system.box.Lx), (self.shrink_time, self.target_length)],
                zero=0
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
