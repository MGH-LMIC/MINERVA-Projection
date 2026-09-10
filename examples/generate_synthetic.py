"""Convenience wrapper; creates only fictitious inputs, never downloads data."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from minerva_projection.synthetic import generate_synthetic

if __name__ == "__main__":
    generate_synthetic(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "synthetic")
