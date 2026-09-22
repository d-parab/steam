"""Command-line interface for STEAM."""

import argparse

from steam.harmonic import harmonic_properties, format_report, write_report


def main():
    parser = argparse.ArgumentParser(
        prog="steam",
        description="Surface Thermodynamics and Entropy Analysis Module",
    )
    parser.add_argument("freq_file", help="text file of frequencies in cm^-1")
    parser.add_argument("-T", "--temperature", type=float, default=298.15,
                        help="temperature in K (default: 298.15)")
    parser.add_argument("-o", "--output", default=None,
                        help="also write results to this file")
    args = parser.parse_args()

    results = harmonic_properties(args.freq_file, args.temperature)
    print(format_report(results, file_path=args.freq_file))

    if args.output:
        write_report(results, args.output, file_path=args.freq_file)
        print(f"\nWritten to {args.output}")


if __name__ == "__main__":
    main()