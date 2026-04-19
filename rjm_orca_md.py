#!/usr/bin/env python3

import numpy as np
import argparse
import os

# get our own functions module
import core_funcs as cf

'''
Simple script for setting up MD in orca

takes a molecule as .xyz, builds an initial configuration

builds an orca input file that should run (lets test that)
'''


# simplified orca input writier (see xyz2orca.py)
def write_orca_input(args, xyz_file):

    base = os.path.splitext(os.path.basename(args.i))[0]
    inp_file = f"{base}_md.inp"

    with open(inp_file, "w") as f:

        f.write(f"! MD {args.method}\n")
        f.write("%pal nprocs {}\nend\n".format(args.nprocs))

        f.write("%md\n")
        f.write(f"  initvel {args.temp}_K\n")
        f.write(f"  timestep {args.timestep}_fs\n")
        f.write(f"  thermostat {args.thermostat} {args.temp}_K timecon {args.timecon}_fs\n")
        f.write(f'  dump position stride {args.dumpfreq} filename "{args.traj}"\n')
        f.write(f"  run {args.nsteps}\n")
        f.write("end\n\n")

        f.write(f"* xyz 0 1\n")
        with open(xyz_file) as xyz:
            lines = xyz.readlines()[2:]
            for l in lines:
                f.write(l)
        f.write("*\n")


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
    parser.add_argument("-ts", "--timestep", type=float, default=0.5)
    parser.add_argument("-tm", "--thermostat", default="berendsen")
    parser.add_argument("-tc", "--timecon", type=float, default=10.0)
    parser.add_argument("-to", "--traj", default="trajectory.xyz")
    parser.add_argument("-df", "--dumpfreq", type=int, default=25)
    parser.add_argument("-ns", "--nsteps", type=int, default=2000)

    parser.add_argument("--nprocs", type=int, default=1)

    args = parser.parse_args()

    atoms, coords = cf.read_xyz(args.i)

    if args.mode == "box":
        atoms_out, coords_out = cf.build_box(atoms, coords, args.n, args.box)
    else:
        atoms_out, coords_out = cf.build_cluster(atoms, coords, args.n, args.dist)

    cf.write_xyz(args.output, atoms_out, coords_out)
    write_orca_input(args, args.output)

if __name__ == "__main__":
    main()