#!/usr/bin/env python3

import numpy as np
import argparse
import os

# get our own functions module
import core_funcs as cf

'''
Simple script for setting up MD in xrb

Note! This is basically the same as the rjm_orca_md one, but running it
directly in XTB seems to make it a little bit faster. 

takes a molecule as .xyz, builds an initial configuration

builds an xtb input file that should run (lets test that)
'''


def write_xtb_input(args, xyz_file):

    base = os.path.splitext(os.path.basename(args.i))[0]
    md_file = f"{base}_md.inp"
    run_file = f"{base}_md.sh"

    # convert steps → time (ps)
    time_ps = args.nsteps * args.timestep / 1000.0

    with open(md_file, "w") as f:

        f.write("$md\n")
        f.write(f"  temp={args.temp}\n")
        f.write(f"  time={time_ps:.3f}\n")
        f.write(f"  step={args.timestep}\n")
        f.write(f"  dump={args.dumpfreq}\n")
        f.write(f"  nvt=true\n")
        f.write(f"  hmass=4\n")
        f.write(f"  shake={args.shake}\n")
        f.write(f"  sccacc={args.sccacc}\n")

        if args.velo:
            f.write(f"  velo=true\n")

        f.write("$end\n")

    with open(run_file, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"xtb {xyz_file} --md --input {md_file} > {base}_md.out\n")

    os.chmod(run_file, 0o755)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", required=True, help="input xyz")
    parser.add_argument("-n", type=int, default=2)
    parser.add_argument("-o", "--output",  default="pre_trajectory.xyz")

    parser.add_argument("--mode", choices=["box","cluster"], default="cluster")
    parser.add_argument("--box", type=float, default=20.0)
    parser.add_argument("--dist", type=float, default=5.0)

    # MD settings
    parser.add_argument("-m", "--method", default="GFN2-XTB")
    parser.add_argument("-t", "--temp", type=float, default=300)
    parser.add_argument("-ts", "--timestep", type=float, default=1.0)
    parser.add_argument("-tm", "--thermostat", default="berendsen")
    parser.add_argument("-tc", "--timecon", type=float, default=10.0)
    parser.add_argument("-to", "--traj", default="trajectory.xyz")
    parser.add_argument("-df", "--dumpfreq", type=int, default=50)
    parser.add_argument("-ns", "--nsteps", type=int, default=2000)

    parser.add_argument("--shake", type=int, default=2)
    parser.add_argument("--sccacc", type=float, default=2.0)
    parser.add_argument("--velo", action="store_true")

    parser.add_argument("--nprocs", type=int, default=1)

    args = parser.parse_args()

    atoms, coords = cf.read_xyz(args.i)

    if args.mode == "box":
        atoms_out, coords_out = cf.build_box(atoms, coords, args.n, args.box)
    else:
        atoms_out, coords_out = cf.build_cluster(atoms, coords, args.n, args.dist)

    cf.write_xyz(args.output, atoms_out, coords_out)
    write_xtb_input(args, args.output)

if __name__ == "__main__":
    main()