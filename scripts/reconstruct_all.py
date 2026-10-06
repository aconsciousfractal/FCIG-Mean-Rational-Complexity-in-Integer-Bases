"""Reconstruct and distinctly check every computational proof premise.

python -B scripts/reconstruct_all.py --work full_normal
python -B scripts/reconstruct_all.py --work full_optimized --optimized
Work names must be fresh. Partial --stage commands do not report completion.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "companion/programs"))
from reconstruction_chain import main

if __name__ == "__main__":
    main()
