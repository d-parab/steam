"""Command-line interface for STEAM."""

import argparse

from steam.harmonic import (
    harmonic_properties,
    format_report as format_harmonic,
    write_report as write_harmonic,
)
from steam.hindered_translation import (
    hindered_translation_properties,
    format_report as format_hindtrans,
    write_report as write_hindtrans,
)
from steam.hindered_rotation import (
    hindered_rotation_properties,
    format_report as format_hindrot,
    write_report as write_hindrot,
)


def main():
    parser = argparse.ArgumentParser(
        prog="steam",
        description="Surface Thermodynamics and Entropy Analysis Module",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # harmonic
    h = sub.add_parser("harmonic",
                       help="harmonic oscillator thermodynamics")
    h.add_argument("freq_file", help="text file of frequencies in cm^-1")
    h.add_argument("-T", "--temperature", type=float, default=298.15,
                   help="temperature in K (default: 298.15)")
    h.add_argument("-o", "--output", default=None,
                   help="also write results to this file")

    # hindered translation
    t = sub.add_parser("hindered-trans",
                       help="hindered translation thermodynamics")
    t.add_argument("contcar", help="VASP CONTCAR file")
    t.add_argument("-W", "--barrier", type=float, required=True,
                   help="translation barrier in eV")
    t.add_argument("-b", "--nn-distance", type=float, required=True,
                   help="nearest neighbour distance in A")
    t.add_argument("-T", "--temperature", type=float, default=298.15,
                   help="temperature in K (default: 298.15)")
    t.add_argument("-o", "--output", default=None,
                   help="also write results to this file")

    # hindered rotation
    r = sub.add_parser("hindered-rot",
                       help="hindered rotation thermodynamics")
    r.add_argument("contcar", help="VASP CONTCAR file")
    r.add_argument("-W", "--barrier", type=float, required=True,
                   help="rotation barrier in eV")
    r.add_argument("-a", "--atom-index", type=int, required=True,
                   help="index of the atom about which the molecule rotates")
    r.add_argument("-n", "--n-minima", type=int, required=True,
                   help="number of equivalent minima in a full rotation")
    r.add_argument("-s", "--symmetry-number", type=int, required=True,
                   help="symmetry number")
    r.add_argument("-T", "--temperature", type=float, default=298.15,
                   help="temperature in K (default: 298.15)")
    r.add_argument("-o", "--output", default=None,
                   help="also write results to this file")

    args = parser.parse_args()

    if args.command == "harmonic":
        results = harmonic_properties(args.freq_file, args.temperature)
        print(format_harmonic(results, file_path=args.freq_file))
        if args.output:
            write_harmonic(results, args.output, file_path=args.freq_file)

    elif args.command == "hindered-trans":
        results = hindered_translation_properties(
            args.barrier, args.contcar, args.nn_distance, args.temperature)
        print(format_hindtrans(results, contcar_file=args.contcar))
        if args.output:
            write_hindtrans(results, args.output, contcar_file=args.contcar)

    elif args.command == "hindered-rot":
        results = hindered_rotation_properties(
            args.barrier, args.contcar, args.atom_index, args.n_minima,
            args.symmetry_number, args.temperature)
        print(format_hindrot(results, contcar_file=args.contcar))
        if args.output:
            write_hindrot(results, args.output, contcar_file=args.contcar)

    if args.output:
        print(f"\nWritten to {args.output}")


if __name__ == "__main__":
    main()