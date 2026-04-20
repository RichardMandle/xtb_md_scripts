#!/usr/bin/env python3

import numpy as np
import mdtraj as md
import core_funcs as cf
import argparse
import matplotlib.pyplot as plt


'''
Analysis pipeline for XTB/Orca MD of small systems

functions:

    part 1 - I/O, trajectory
parse_md_input              - read the xtb/orca md input file and - more or less - workout the timesteps per frame in the trajectory.
parse_single_xyz            - reads the .xyz file containing coordinates of one molecule
get_natoms_from_xyz         - count the number of atoms in one molecule
build_topology              - build an mdtraj compatible topology 
load_trajectory             - load the .xyz trajectory as mdtraj traj object 

    part 2 - basic analysis

compute_com                 - compute centre-of-mass (COM) for each molecule
compute_pair_distnces       - compute pairwise COM-COM distance per frame
compute_min_distances       - compute the distance between molecules on a minimum atom-atom basis
compute_interaction_mask    - mask out interactions beyond a cutoff distance (controleld with -c)
compute_orientation_when... -


    part 3 - angle based interactions
get_directors               - use mdtraj to compute the vector describing the molecules orientation
compute_pair_angles         - calculate angles between all pairs of molecules using the vectors from get_directors
fix_director_signs          - keep the director signs constant; stops the vector description of molecular orientation "flipping"

    part n - plotting tools
plot_distances              - plot pair COM-COM distances per frame
plot_orientations           - plot pair cos(theta) per frame
plot_orientation_hist       - histogram of cos(thera) per pair over all frames

TO DO:
Different analysis tools. COM-COM is too blunt. Would be good to decompose that into x/y/z
in thre frame of the molecule (e.g. mol 0 interacts with mol 1; its displaced along X by...)

'''

# Part 1
def parse_md_input(inp_file):
    step = None
    dump = None

    with open(inp_file) as f:
        for line in f:
            line = line.strip()

            if line.startswith("step"):
                step = float(line.split("=")[1])

            elif line.startswith("dump"):
                dump = int(line.split("=")[1])

    if step is None or dump is None:
        raise ValueError(f"Could not parse step/dump from {inp_file}")

    timestep_fs = step * dump

    print(f"[INFO] MD input:")
    print(f"       step = {step} fs")
    print(f"       dump = {dump}")
    print(f"       frame timestep = {timestep_fs} fs")

    return timestep_fs
    
def parse_single_xyz(xyz_file):
    # reads a single molecule in an xyz file
    atoms, coords = cf.read_xyz(xyz_file)
    nat_per_mol = len(atoms)

    print(f"[INFO] Single molecule: atoms per molecule = {nat_per_mol}")

    return nat_per_mol, atoms, coords

def get_natoms_from_xyz(traj_file):
    # reads an xyz trajectory and gets the number of atoms
    with open(traj_file) as f:
        first_line = f.readline()

    try:
        n_atoms = int(first_line.strip())
    
    except Exception:
        raise ValueError(f"Could not read atom count from {traj_file}")

    print(f"[INFO] Trajectory: atoms per frame = {n_atoms}")

    return n_atoms

def build_topology(atoms_single, nmol):
    # builds a topology object so we can use mdtraj to load the trajectory
    top = md.Topology()
    chain = top.add_chain()

    for i in range(nmol):
        res = top.add_residue(f"MOL{i}", chain)

        for atom in atoms_single:
            element = md.element.get_by_symbol(atom)
            top.add_atom(atom, element, res)

    return top


def load_trajectory(traj_file, atoms_single, nat_per_mol):
    # use the built topology to load the trajectory
    n_atoms = get_natoms_from_xyz(traj_file)

    if n_atoms % nat_per_mol != 0:
        raise ValueError("Trajectory atom count not divisible by molecule size")

    nmol = n_atoms // nat_per_mol

    print(f"[INFO] Derived system: number of molecules = {nmol}")

    top = build_topology(atoms_single, nmol)
    traj = md.load_xyz(traj_file, top=top)

    print(f"\n[INFO] Loaded trajectory:")
    print(f"       Frames = {traj.n_frames}")
    print(f"       Atoms  = {traj.n_atoms}")

    return traj, nmol

# part 2 - basic analysis
def compute_com(traj):
    coms = []

    for res in traj.topology.residues:
        inds = [atom.index for atom in res.atoms]
        coords = traj.xyz[:, inds, :]  
        com = coords.mean(axis=1)      
        coms.append(com)

    return np.array(coms)  # tuple of (nmol, nframes, 3)

def compute_pair_distances(coms):
    nmol, nframes, _ = coms.shape

    pair_dists = {}

    for i in range(nmol):
        for j in range(i+1, nmol):
            rij = np.linalg.norm(coms[i] - coms[j], axis=1)
            pair_dists[(i, j)] = rij  # (nframes,)

    return pair_dists
    
def compute_min_distances(traj):
    """
    compute minimum atom–atom distance between molecule pairs per frame.

    return:
    pair_min_dists: dict[(i,j)] -> (n_frames,)
    """
    pair_min_dists = {}

    mol_atoms = [
        [atom.index for atom in res.atoms]
        for res in traj.topology.residues
    ]

    nmol = len(mol_atoms)
    coords = traj.xyz

    for i in range(nmol):
        for j in range(i+1, nmol):

            inds_i = mol_atoms[i]
            inds_j = mol_atoms[j]

            xyz_i = coords[:, inds_i, :]
            xyz_j = coords[:, inds_j, :]

            d = np.linalg.norm(
                xyz_i[:, :, None, :] - xyz_j[:, None, :, :],
                axis=-1
            )

            d_min = d.min(axis=(1, 2))

            pair_min_dists[(i, j)] = d_min

    return pair_min_dists
    
def compute_interaction_mask(args, pair_min_dists, cutoff=0.5):
    """
    cutoff in nm
    """
    masks = {}

    for k, d in pair_min_dists.items():
        masks[k] = d < cutoff

    return masks
    
def compute_interaction_fraction(masks):
    print("\n[ANALYSIS] Interaction fraction:")

    for k, mask in masks.items():
        frac = mask.sum() / len(mask)
        print(f"{k}: {frac*100:.1f}%")
        
def compute_orientation_when_interacting(pair_angles, masks):
    print("\n[ANALYSIS] Orientation when interacting:")

    for k in pair_angles.keys():
        c = pair_angles[k]
        mask = masks[k]

        if np.any(mask):
            mean = c[mask].mean()
            std = c[mask].std()
            print(f"{k}: mean cosθ = {mean:.2f} ± {std:.2f}")
        else:
            print(f"{k}: no interactions")
            
            
#part 3 - angle analysis
 
def get_directors(traj, indices='residues', care_about_polar = True):
    # return a vector describing the orientation of each molecule in each frame
    directors = md.compute_directors(traj, indices=indices)
    directors = np.transpose(directors, (1, 0, 2))  # (n_mols, n_frames, 3)
    
    if care_about_polar:
        directors = fix_director_signs(directors)
        
        # delete this later RJM TO DO
        print("\n[DEBUG] director continuity check:")
        for m in range(min(2, directors.shape[0])):
            dots = np.sum(directors[m, 1:] * directors[m, :-1], axis=1)
            print(f"mol {m}: min dot(prev, curr) = {dots.min():.3f}")
            
    return directors

def compute_pair_angles(traj, directors):
    n_mols, n_frames, _ = directors.shape

    pair_angles = {}

    for i in range(n_mols):
        for j in range(i+1, n_mols):
            cos_theta = np.sum(directors[i, :, :] * directors[j, :, :], axis=1)
            pair_angles[(i, j)] = cos_theta

    return pair_angles
    
def fix_director_signs(directors):
    """
    If we care about polarity then we need to have time-continuity of director vectors.
    directors: (n_mols, n_frames, 3)
    """
    n_mols, n_frames, _ = directors.shape

    for m in range(n_mols):
        for t in range(1, n_frames):
            if np.dot(directors[m, t], directors[m, t-1]) < 0:
                directors[m, t] *= -1

    return directors

# part n - plotting
def plot_distances(pair_dists, timestep_fs=1.0, outfile="distances.png"):
    '''
    very basic plotter for compute_pair_distances data
    '''
    plt.figure()

    for (i, j), d in pair_dists.items():
        t = np.arange(len(d)) * timestep_fs / 1000.0  # ps
        plt.plot(t, d, label=f"{i}-{j}")

    plt.xlabel("Time (ps)")
    plt.ylabel("COM distance (nm)")
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(outfile, dpi=300)
    plt.show()

def plot_orientations(pair_angles, timestep_fs=1.0, outfile="orientations.png", ignore_polar = False):
    plt.figure()

    for (i, j), c in pair_angles.items():
        t = np.arange(len(c)) * timestep_fs / 1000.0  # ps
        if ignore_polar:
            c = np.abs(c)
        plt.plot(t, c, label=f"{i}-{j}", alpha=0.8)

    plt.xlabel("Time (ps)")
    plt.ylabel("cos(theta)")
    plt.ylim(-1.1, 1.1)
    plt.legend()
    plt.tight_layout()

    plt.savefig(outfile, dpi=300)
    plt.show()
    
def plot_orientation_hist(pair_angles, masks, outfile="orientation_hist.png",):
    plt.figure()

    for k in pair_angles.keys():
        c = pair_angles[k]
        mask = masks[k]

        if np.any(mask):
            plt.hist(
                c[mask],
                bins=30,
                alpha=0.5,
                label=f"{k}"
            )

    plt.xlabel("cos(theta)")
    plt.ylabel("count")
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(outfile, dpi=300)
    plt.show()
    
    
def main():
    parser = argparse.ArgumentParser(description="Process an XTB/ORCA trajectory")

    parser.add_argument("-i", "--input_xyz", required=True, help="Single molecule XYZ (template)")
    parser.add_argument("-t", "--traj", default="xtb_trajectory.xyz", help="Trajectory XYZ file")
    parser.add_argument("-m", "--md_input", required=True, help="xTB/ORCA MD input file (.inp)" )
    
    # turn on different analysis routines  
    parser.add_argument("-a", "--angle", action="store_true", help="Run orientation (angle) analysis")
    parser.add_argument("-d", "--distance", action="store_true", help="Run distance analysis")
    parser.add_argument("-ad", "--interaction", action="store_true", help="Run interaction vs angle analysis")

    # tunable parameters stuff  
    parser.add_argument("-c", "--int_cutoff", default=0.5, help="Minimum atom-atom distance (in nm) required for interaction")
    
    args = parser.parse_args()

    timestep_fs = parse_md_input(args.md_input)
    
    nat_per_mol, atoms_single, coords_single = parse_single_xyz(args.input_xyz)
    traj, nmol = load_trajectory(args.traj, atoms_single, nat_per_mol)
        
    if args.distance: # Distance stuff    
        coms = compute_com(traj)
        pair_dists = compute_pair_distances(coms)
        plot_distances(pair_dists, timestep_fs=timestep_fs)
    
    if args.angle: # Orientation stuff
        directors = get_directors(traj)
        pair_angles = compute_pair_angles(traj, directors)
        plot_orientations(pair_angles, timestep_fs=timestep_fs)
    
    if args.interaction: # Cos(theta) vs interaction histogram
        if 'pair_angles' not in locals():
            directors = get_directors(traj)
            pair_angles = compute_pair_angles(traj, directors)
            
        pair_min_dists = compute_min_distances(traj)
        masks = compute_interaction_mask(args, pair_min_dists, cutoff=args.int_cutoff)
        compute_interaction_fraction(masks)
        plot_orientation_hist(pair_angles, masks)
        
        
if __name__ == "__main__":
    main()
