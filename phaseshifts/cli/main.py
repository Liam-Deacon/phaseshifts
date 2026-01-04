import argparse
import sys
from phaseshifts import phsh


def main():
    parser = argparse.ArgumentParser(
        prog="phaseshifts",
        description="Phase shifts and atomic orbital library generator CLI.",
    )
    subparsers = parser.add_subparsers(dest="command", required=False)

    # phsh subcommand (default)
    phsh_parser = subparsers.add_parser(
        "phsh",
        help="Generate phase shifts (default)",
        add_help=False,  # We'll let phsh.py handle its own help
    )
    phsh_parser.add_argument("phsh_args", nargs=argparse.REMAINDER)

    # atorb subcommand
    atorb_parser = subparsers.add_parser("atorb", help="Generate atomic orbital input files for all known elements")
    atorb_parser.add_argument(
        "--output-dir",
        type=str,
        default="./atorb_lib",
        help="Directory to write all generated atorb input files",
    )
    atorb_parser.add_argument(
        "--rel",
        action="store_true",
        default=False,
        help="Use relativistic calculation (default: False)",
    )
    atorb_parser.add_argument("--ngrid", type=int, default=1000, help="Number of grid points (default: 1000)")
    atorb_parser.add_argument(
        "--method",
        type=str,
        default="HF",
        help="Exchange-correlation method (default: HF)",
    )

    args, unknown = parser.parse_known_args()

    # Default to phsh if no subcommand is given
    if args.command is None or args.command == "phsh":
        # Forward all args to phsh.main
        sys.argv = ["phsh.py"] + (getattr(args, "phsh_args", []) or []) + unknown
        return phsh.main()

    elif args.command == "atorb":
        import os

        os.makedirs(args.output_dir, exist_ok=True)
        from phaseshifts.atorb import elements_dict, Atorb

        generated = []
        for symbol in elements_dict:
            try:
                filename = Atorb.gen_input(
                    symbol,
                    rel=args.rel,
                    ngrid=args.ngrid,
                    method=args.method,
                    filename=os.path.join(args.output_dir, f"atorb_{symbol}.txt"),
                    output=os.path.join(args.output_dir, f"at_{symbol}.i"),
                )
                print(f"Generated: {filename}")
                generated.append(filename)
            except Exception as e:
                print(f"Failed to generate for {symbol}: {e}", file=sys.stderr)
        print(f"\nGenerated {len(generated)} atorb input files in {args.output_dir}")
        return 0


if __name__ == "__main__":
    sys.exit(main())
