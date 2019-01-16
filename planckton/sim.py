import os
import logging
import hoomd.deprecated
import hoomd.md
import hoomd.dump
import hoomd.data
from cme_utils.manip.convert_rigid import init_wrapper
from cme_utils.manip.ff_from_foyer import set_coeffs


class Simulation:
    def __init__(
        self,
        input_xml,
        kT,
        e_factor=1.0,
        tau=5.0,
        gsd_write=1e6,
        log_write=1e5,
        shrink_time=1e6,
        shrink_factor=5,
        n_steps=1e3,
        mode="gpu",
    ):
        self.input_xml = input_xml
        self.e_factor = e_factor
        self.tau = tau
        self.kT = kT
        self.gsd_write = gsd_write
        self.log_write = log_write
        self.shrink_time = shrink_time
        self.shrink_factor = shrink_factor
        self.n_steps = n_steps
        self.mode = mode

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
            integrator_mode = hoomd.md.integrate.mode_standard(dt=0.0001)
            rigid = hoomd.group.rigid_center()
            nonrigid = hoomd.group.nonrigid()
            both_group = hoomd.group.union("both", rigid, nonrigid)
            all_particles = hoomd.group.all()
            integrator = hoomd.md.integrate.nvt(
                group=both_group, tau=self.tau, kT=self.kT
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
            desired_box_dim = system.box.Lx / self.shrink_factor
            size_variant = hoomd.variant.linear_interp(
                [(0, system.box.Lx), (self.shrink_time, desired_box_dim)], zero=0
            )
            hoomd.update.box_resize(L=size_variant)
            hoomd.run_upto(self.shrink_time)
            # After shrinking, reset velocities
            integrator.randomize_velocities(seed=42)
            integrator_mode.set_params(dt=0.0001)
            try:
                hoomd.run_upto(self.n_steps+1, limit_multiple=self.gsd_write)
            except hoomd.WalltimeLimitReached:
                pass
            finally:
                gsd_restart.write_restart()


if __name__ == "__main__":
    my_sim = Simulation("init.hoomdxml", kT=3.0, gsd_write=1e5, log_write=1e5)
    my_sim.run()
