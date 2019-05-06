import os
import simulation
import numpy as np
import pandas as pd
from simulation.parsers.parse_info import parse_info
from simulation.units import gadget_time_units
from simulation.util import get_sim_name, make_info_monotonic_again
from simulation.derotate_simulation import rotate_vec
from astropy.table import Table
import matplotlib.pyplot as plt
from pyquaternion import Quaternion
import tqdm
import quaternion
from simulation.NFW2acc import nfw
import argparse


def Omega(v, a):
    return np.cross(v, a)/((np.linalg.norm(v, axis=1)**2)[:, np.newaxis])


def get_acceleration(trace, nfw=nfw):
    pos = np.vstack([trace.x, trace.y, trace.z]).T
    acc_rec = nfw.acc(pos)
    adhoc = np.vstack([trace.ax_adhoc, trace.ay_adhoc, trace.az_adhoc]).T
    real_adhoc = adhoc * np.linalg.norm(acc_rec, axis=1)[:, np.newaxis]
    tot_acc = acc_rec + real_adhoc
    return tot_acc


def get_Omega(trace, nfw=nfw):
    acc = get_acceleration(trace, nfw)
    v_all = np.vstack([trace.vx, trace.vy, trace.vz]).T
    return Omega(v_all, acc)


def get_quaternions(trace):
    quat = np.vstack([trace.quat_w, trace.quat_x, trace.quat_y, trace.quat_z]).T
    return quaternion.as_quat_array(quat)


def evolve_quaternion_q(omega, timesteps, q0=np.quaternion(1,0,0,0)):
    """Evolve quaternion. Return quaternion"""
    q_list = list()
    q1 = q0
    for i, dt in enumerate(tqdm.tqdm(timesteps[1:])):
        quat_deriv = 0.5 * np.quaternion(0, *omega[i]) * q1
        q1 += quat_deriv * dt
        q_list.append(q1/q1.norm())
    return np.squeeze(np.vstack(q_list))


def compute_quaternion(trace, dt):
    """Return quaternion"""
    omega = get_Omega(trace)
    q = evolve_quaternion_q(omega, dt)
    return q


def compute_and_write_quaternion(sim_path, full_sim=False, force_recovery=False):
    sim = simulation.Simulation(sim_path)
    # print(sim.trace.head())
    if 'quat_w' in sim.trace:
        if force_recovery:
            print("Recomputing quaternion even if I already have it in the trace file...")
            q = compute_quaternion(sim.trace, sim.trace.dt)

        else:
            print("Quaternion already there, just saving it...")
            q = quaternion.as_quat_array(np.vstack([sim.trace.quat_w, sim.trace.quat_x, sim.trace.quat_y, sim.trace.quat_z]).T)
    else:
        print("Recovering quaternion...")
        # I need dt which is in info.txt
        info = parse_info(os.path.join(sim_path, "info.txt"))
        if not info.step.is_monotonic_increasing:
            print("Adjusting info...")
            info_orig = info.copy()
            info = make_info_monotonic_again(info_orig)
        # print(info.head())
        assert len(info) == len(sim.trace), (len(info), len(sim.trace))

        q = compute_quaternion(sim.trace, info.dt)

    if full_sim:
        outname = get_sim_name(sim_path) + "_quat_fullsim.fits"
        location = slice(None) # TODO check
    else:
        # Skipping the first row because the first time is duplicated? FIXME
        locations = np.digitize(sim.times_header, sim.trace.t[1:], right=True)
        q = q[locations]
        outname = get_sim_name(sim_path) + "_quat.fits"

    q_code = q.copy()

    # Add quaternion of the velocity
    vp0 = np.array([sim.trace.vx[0], sim.trace.vy[0], sim.trace.vz[0]])
    print("Vp = {}".format(vp0))
    # This is the bisection formula applied to the angle between V_p and -y
    cosv0 = -vp0[1]/np.linalg.norm(vp0)
    quat_vp0 = np.quaternion(-np.sqrt((1+cosv0)/2), 0, 0, np.sqrt((1-cosv0)/2))
    print("Adding initial quaternion (V_p -> -y)...")
    q *= quat_vp0

    # Get Omega
    if 'omega_x' not in sim.trace or force_recovery:
        Omega = get_Omega(sim.trace)
    else:
        Omega = np.hstack([sim.trace[['omega_x', 'omega_y', 'omega_z']]])
    Omega = Omega[locations]

    # Get omega_mb
    if 'omega_mb_x' not in sim.trace or force_recovery:
        omega_mb = rotate_vec(Omega, q_code.conj())  # This is what is done in the code. I doubt it is correct
    else:
        omega_mb = np.hstack([sim.trace[['omega_mb_x', 'omega_mb_y', 'omega_mb_z']]])[locations]
        # omega_mb = rotate_vec(omega_mb, quat_vp0)

    Omega = rotate_vec(Omega, quat_vp0)
    omega_mb = rotate_vec(omega_mb, quat_vp0)


    # write quaternion
    tbl = Table(np.hstack([quaternion.as_float_array(q), Omega, omega_mb]),
                names=["q_w", "q_x", "q_y", "q_z", 'omega_x', 'omega_y', 'omega_z', 'omega_mb_x', 'omega_mb_y', 'omega_mb_z'],
                meta={'V_P0X': sim.trace.vx[0], "V_P0X": sim.trace.vy[0], "V_P0Z": sim.trace.vz[0]})

    tbl.write(outname, overwrite=True)
    return sim


def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("sim_path", help="Simulation path")
    parser.add_argument("--full-sim", help="Save quaternion for all the timesteps of the simulation", action="store_true")
    parser.add_argument("--force-recovery", "-f", help="Recover quaternion even if it is already in the trace file", action="store_true")
    # TODO
    parser.add_argument("--outname", "-o", help="Output file name", default=None)
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
    compute_and_write_quaternion(args.sim_path, args.full_sim, args.force_recovery)


if __name__ == '__main__':
    main()


