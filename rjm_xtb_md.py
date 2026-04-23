#!/usr/bin/env python3

import json
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

def resolve_base(args):
    '''
    NEW
    properly resolve basenames for when we are using replicate simulations
    '''
    if args.prefix:
        return args.prefix
    return os.path.splitext(os.path.basename(args.output))[0]
    
def ensure_xyz_suffix(filename):
    '''
    NEW
    Seems XTB might error if we try to pass an xyz file without the .xyz extension
    this can happen if the -o flag doesn't end .xyz; so lets enforce it.
    '''
    if not filename.lower().endswith(".xyz"):
        filename += ".xyz"
    return filename
    
def write_xtb_input(args, xyz_file, base):

    md_file = f"{base}_md.inp"
    run_file = f"{base}_md.sh"

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

    method_flag = ""
    if args.method.upper() == "GFN2-XTB":
        method_flag = "--gfn 2"
    elif args.method.upper() == "GFN1-XTB":
        method_flag = "--gfn 1"
        
# For local usage we write a .sh file we can execute with zsh 
    if args.enviro.lower() == 'local':
        with open(run_file, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("ulimit -s unlimited\n")
            f.write(f"export OMP_NUM_THREADS={args.nprocs},1\n")
            f.write(f"export MKL_NUM_THREADS={args.nprocs}\n")
            f.write("export OMP_MAX_ACTIVE_LEVELS=1\n")
            f.write("export OMP_STACKSIZE=4G\n")
            f.write(f"xtb {xyz_file} {method_flag} -P {args.nprocs} --md --input {md_file} > {base}_md.out\n")
            
# For HPC usage assume we'll send the .sh file to slurm via sbatch; so write the script appropriately.
# Use $TMP_LOCAL for on-node NVME scatch during run.
    if args.enviro.lower() == 'hpc':
        with open(run_file, "w") as f:
            f.write(f"""#!/bin/bash
#SBATCH --job-name={base}_xtb_md
#SBATCH --output={base}_md.out
#SBATCH --error={base}.err
#SBATCH --ntasks={args.nprocs}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem={args.nprocs*4}G
#SBATCH --time=2-00:00:00

module add openmpi
module load xtb

working_dir=$(pwd)

cp -a ./* $TMP_LOCAL/
cd $TMP_LOCAL

xtb {xyz_file} {method_flag} -P {args.nprocs} --md --input {md_file} > {base}_md.out

cp -a ./* "$working_dir"/
cd "$working_dir"

""")
        os.chmod(run_file, 0o755)


def write_manifest(args, base):
    '''
    NEW
    Log everything in a json manifest
    '''
    manifest = vars(args).copy()
    with open(f"{base}_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", required=True, help="input xyz")
    parser.add_argument("-n", type=int, default=2)
    parser.add_argument("-o", "--output", default="pre_trajectory.xyz", help="name of the coordinate file for the simulation initial configuration (should end .xyz!)")
    parser.add_argument("-prefix", default=None, help="base name for generated files")
    parser.add_argument("-seed", type=int, default=None, help="RNG seed for reproducible packing")
    parser.add_argument("-rep", type=int, default=0, help="replicate index")

    parser.add_argument("-mode", choices=["box", "cluster"], default="cluster")
    parser.add_argument("-enviro", choices=["local", "hpc"], default="local")
    parser.add_argument("-box", type=float, default=50.0)
    parser.add_argument("-dist", type=float, default=5.0)

    parser.add_argument("-m", "--method", default="GFN2-XTB")
    parser.add_argument("-t", "--temp", type=float, default=300.0)
    parser.add_argument("-ts", "--timestep", type=float, default=1.0)
    parser.add_argument("-df", "--dumpfreq", type=int, default=50)
    parser.add_argument("-ns", "--nsteps", type=int, default=2000)

    parser.add_argument("-shake", type=int, default=2)
    parser.add_argument("-sccacc", type=float, default=2.0)
    parser.add_argument("-velo", action="store_false")
    parser.add_argument("-nprocs", type=int, default=4)

    args = parser.parse_args()
    args.output = ensure_xyz_suffix(args.output) 
    
    if not os.path.isfile(args.i):
        raise FileNotFoundError(f"Input file not found: {args.i}")
    if args.n < 1:
        raise ValueError("n must be >= 1")
    if args.timestep <= 0:
        raise ValueError("timestep must be > 0")
    if args.nsteps <= 0:
        raise ValueError("nsteps must be > 0")
    if args.dumpfreq <= 0:
        raise ValueError("dumpfreq must be > 0")

    seed = None if args.seed is None else int(args.seed) + int(args.rep)
    rng = None if seed is None else __import__("numpy").random.default_rng(seed)

    atoms, coords = cf.read_xyz(args.i)

    if args.mode == "box":
        atoms_out, coords_out = cf.build_box(atoms, coords, args.n, args.box, rng=rng)
    else:
        atoms_out, coords_out = cf.build_cluster(atoms, coords, args.n, args.dist, rng=rng)

    base = resolve_base(args)

    outdir = os.path.dirname(args.output)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    cf.write_xyz(args.output, atoms_out, coords_out)
    write_xtb_input(args, args.output, base)
    write_manifest(args, base)


if __name__ == "__main__":
    main()