#!/usr/bin/env python3

import numpy as np

'''
core_funcs - just some core functions we'll use across different scripts.

list:
read_xyz        -   basic .xyz file reader
write_xyz       -   basic .xyz file writer
random_rotation -   randomly rotate molecules
centre          -   literally just centre the coordinates  
'''


# Functions for handling .xyz files
def read_xyz(filename):
    with open(filename) as f:
        lines = f.readlines()
    n = int(lines[0])
    atoms, coords = [], []
    for line in lines[2:2+n]:
        a, x, y, z = line.split()
        atoms.append(a)
        coords.append([float(x), float(y), float(z)])
    return atoms, np.array(coords)
    

def write_xyz(filename, atoms, coords):
    with open(filename, "w") as f:
        f.write(f"{len(atoms)}\n\n")
        for a, (x,y,z) in zip(atoms, coords):
            f.write(f"{a} {x:.6f} {y:.6f} {z:.6f}\n")


# Functions for simulation box packing, molecule manipuation etc.
def random_rotation():
    u1, u2, u3 = np.random.rand(3)
    q = [
        np.sqrt(1-u1)*np.sin(2*np.pi*u2),
        np.sqrt(1-u1)*np.cos(2*np.pi*u2),
        np.sqrt(u1)*np.sin(2*np.pi*u3),
        np.sqrt(u1)*np.cos(2*np.pi*u3)
    ]
    w,x,y,z = q
    return np.array([
        [1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
        [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)],
        [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)]
    ])

def center(coords):
    return coords - coords.mean(axis=0)
    


# We need a whole bunch of functions for packing into the md box
# atomic clash cutoff (Å) — deliberately conservative
CLASH_CUTOFF = 2.5

def atoms_clash(existing, new, cutoff=CLASH_CUTOFF):
    if len(existing) == 0:
        return False
    d = np.linalg.norm(existing[:,None,:] - new[None,:,:], axis=-1)
    return np.any(d < cutoff)

def estimate_span(coords):
    """max molecular extent (used for box padding)"""
    return np.max(coords.max(axis=0) - coords.min(axis=0))

def build_box(atoms, coords, n, box=20.0):
    coords = center(coords)

    all_coords = []
    all_atoms = []

    span = estimate_span(coords)
    margin = span + 2.0

    placed = 0
    max_attempts = 2000
    attempts = 0

    while placed < n:

        R = random_rotation()
        t = np.random.uniform(margin, box - margin, 3)

        new = coords @ R.T + t

        existing = np.vstack(all_coords) if all_coords else []

        if not atoms_clash(existing, new):
            all_coords.append(new)
            all_atoms.extend(atoms)
            placed += 1
            attempts = 0

        else:
            attempts += 1

        if attempts > max_attempts:
            raise RuntimeError(
                f"Failed to place molecule {placed+1} — increase box size"
            )

    return all_atoms, np.vstack(all_coords)

# We can specify a distance so molecules aren't too far appart and so spend a lot of
# steps finding each other (or never doing so)
def build_cluster(atoms, coords, n, dist=5.0, max_attempts = 2000):
    coords = center(coords)

    all_atoms = []
    all_coords = []

    all_atoms.extend(atoms)
    all_coords.append(coords)

    placed = 1 
    attempts = 0

    while placed < n:

        R = random_rotation()
        direction = np.random.randn(3)
        direction /= np.linalg.norm(direction)

        new = coords @ R.T + direction * dist

        existing = np.vstack(all_coords)

        if not atoms_clash(existing, new):
            all_atoms.extend(atoms)
            all_coords.append(new)
            placed += 1
            attempts = 0  # reset after success

        else:
            attempts += 1

        # if it struggles, keep trying with looser bounds.
        if attempts > max_attempts:
            dist *= 1.2
            attempts = 0

    return all_atoms, np.vstack(all_coords)
