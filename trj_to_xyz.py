#!/usr/bin/env python3

import argparse

'''
trj_to_xyz.py

little script that converts an xtb trajectory .trj to xyz so we can easily do a 
visual check on it.

For actual analysis, use the trj file directly.
'''

def read_xyz_atoms(xyz_file):
    with open(xyz_file) as f:
        lines = [line.rstrip() for line in f]

    try:
        nat = int(lines[0].strip())
    except Exception as e:
        raise ValueError(f"Could not read atom count from {xyz_file}") from e

    atoms = []
    for line in lines[2:2 + nat]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Malformed XYZ line in {xyz_file}: {line}")
        atoms.append(parts[0])

    if len(atoms) != nat:
        raise ValueError(f"Expected {nat} atoms in {xyz_file}, found {len(atoms)}")

    return atoms


def looks_like_coord_line(line):
    parts = line.split()
    if len(parts) < 4:
        return False
    try:
        float(parts[-3])
        float(parts[-2])
        float(parts[-1])
        return True
    except ValueError:
        return False


def parse_xtb_trj(trj_file, nat):
    with open(trj_file) as f:
        lines = [line.rstrip("\n") for line in f]

    frames = []
    current = []

    for line in lines:
        stripped = line.strip()

        if not stripped:
            continue

        if looks_like_coord_line(stripped):
            parts = stripped.split()
            x, y, z = map(float, parts[-3:])
            current.append((x, y, z))

            if len(current) == nat:
                frames.append(current)
                current = []
        else:
            # ignore headers / energies / metadata lines
            continue

    if current:
        raise ValueError(
            f"Found incomplete final frame in {trj_file}: "
            f"{len(current)} coordinates, expected {nat}"
        )

    if not frames:
        raise ValueError(f"No coordinate frames found in {trj_file}")

    return frames


def write_multixyz(out_file, atoms, frames):
    nat = len(atoms)
    with open(out_file, "w") as f:
        for i, frame in enumerate(frames):
            f.write(f"{nat}\n")
            f.write(f"Frame {i + 1}\n")
            for atom, (x, y, z) in zip(atoms, frame):
                f.write(f"{atom:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert xTB xtb.trj to a multi-frame XYZ using atom labels from an input XYZ."
    )
    parser.add_argument("-i", "--input_xyz", required=True, help="Input XYZ used for xTB MD")
    parser.add_argument("-t", "--traj", default="xtb.trj", help="xTB trajectory file")
    parser.add_argument("-o", "--output", default="xtb_trajectory.xyz", help="Output multi-frame XYZ")
    args = parser.parse_args()

    atoms = read_xyz_atoms(args.input_xyz)
    frames = parse_xtb_trj(args.traj, len(atoms))
    write_multixyz(args.output, atoms, frames)

    print(f"Wrote {len(frames)} frames to {args.output}")


if __name__ == "__main__":
    main()