"""Command-line entry point; all file access is caller-directed and local."""
import argparse
import json
from pathlib import Path
from .pipeline import read_csv, run_analysis, write_outputs
from .synthetic import generate_synthetic


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Project supplied literature directions onto cohort data")
    commands = parser.add_subparsers(dest="command", required=True)
    synthetic = commands.add_parser("synthetic", help="Generate clearly fictitious input CSVs")
    synthetic.add_argument("--output", type=Path, required=True)
    synthetic.add_argument("--seed", type=int, default=17)
    run = commands.add_parser("run", help="Analyze local user-supplied CSVs")
    for name in ["participants", "abundance", "votes", "output"]:
        run.add_argument(f"--{name}", type=Path, required=True)
    run.add_argument("--exclude-publications", type=Path)
    run.add_argument("--min-group-size", type=int, default=10)
    run.add_argument("--min-reference-size", type=int, default=100)
    run.add_argument("--display-genera", type=int, default=28)
    run.add_argument("--min-landscape-sources", type=int, default=2)
    args = parser.parse_args(argv)
    try:
        if args.command == "synthetic":
            generate_synthetic(args.output, args.seed)
            print("Generated entirely synthetic input files.")
            return 0
        excluded = set()
        if args.exclude_publications:
            x = read_csv(args.exclude_publications)
            if "publication_id" not in x:
                raise ValueError("Exclusion CSV requires publication_id")
            excluded = set(x.publication_id)
        tables, summary = run_analysis(read_csv(args.participants), read_csv(args.abundance),
            read_csv(args.votes), excluded_publications=excluded,
            min_group_size=args.min_group_size, min_reference_size=args.min_reference_size,
            display_genera=args.display_genera, min_landscape_sources=args.min_landscape_sources)
        write_outputs(args.output, tables, summary)
        print(json.dumps({"status": summary["status"], "eligible_tasks": summary["eligible_tasks"]}))
        return 0 if summary["status"] == "complete" else 2
    except (ValueError, OSError) as error:
        parser.exit(2, f"Input/output error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
