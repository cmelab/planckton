import logging
import os

import hoomd.data
import hoomd.dump
import hoomd.md
from mbuild.formats.hoomd_simulation import create_hoomd_simulation

from planckton.utils.rigid import init_rigid
from planckton.utils.utils import set_coeffs


class Simulation:
    def __init__(
        self,
        typed_system,
        kT,
        e_factor=1.0,
        rigid_inds=None,
        rigid_typeids=None,
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
        self.system = typed_system
        self.e_factor = e_factor
        self.rigid_inds = rigid_inds
        self.rigid_typeids = rigid_typeids
        self.tau = tau
        self.kT = kT
        self.gsd_write = gsd_write
        self.log_write = log_write
        self.shrink_time = shrink_time
        self.shrink_factor = shrink_factor  # this isn't used anywhere??
        self.shrink_kT_reduced = shrink_kT_reduced
        self.n_steps = n_steps
        self.dt = dt
        self.mode = mode
        self.target_length = target_length

    def run(self):
        hoomd_args = f"--single-mpi --mode={self.mode}"
        sim = hoomd.context.initialize(hoomd_args)
        if self.rigid_system is not None:
            both, snap, sim, ref_values = init_rigid(
                self.rigid_inds, self.rigid_typeids, self.system, sim
            )
        else:
            with sim:
                hoomd.util.quiet_status()
                hoomd_objects, ref_values = create_hoomd_simulation(
                    self.system, auto_scale=True
                )
                snap = hoomd_objects[0]
                hoomd.util.unquiet_status()

        with sim:
            if self.target_length is not None:
                self.target_length /= ref_values.distance

            if self.e_factor != 1:
                logging.info("Scaling LJ coeffs by e_factor")
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
                self.target_length = snap.box.Lx  # should be /scale_factor?
            size_variant = hoomd.variant.linear_interp(
                [(0, snap.box.Lx), (self.shrink_time, self.target_length)],
                zero=0,
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
