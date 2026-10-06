"""Build the standalone paper with an existing pdflatex installation."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
JOB = "Mean-Rational-Complexity-in-Integer-Bases"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    paper = ROOT/"paper"
    out = (args.output_dir or paper/"build").resolve()
    print(out, flush=True)
    out.mkdir(parents=True, exist_ok=True)
    exe = os.environ.get("PDFLATEX") or shutil.which("pdflatex")
    if not exe:
        raise SystemExit("PDF_BUILD_FAIL: existing pdflatex required; set PDFLATEX or PATH")
    version = subprocess.run([exe, "--version"], capture_output=True, text=True, check=True)
    flags = ["-no-shell-escape", "-interaction=nonstopmode", "-halt-on-error"]
    if "miktex" in (version.stdout + version.stderr).lower():
        flags.insert(0, "--disable-installer")
    command = [exe, *flags,
               "-output-directory=" + str(out), JOB + ".tex"]
    for _ in range(2):
        subprocess.run(command, cwd=paper, check=True)
    produced = out/(JOB + ".pdf")
    if not produced.is_file():
        raise SystemExit("PDF_BUILD_FAIL: expected PDF missing")
    dest = paper/(JOB + ".pdf")
    if produced.resolve() != dest.resolve():
        shutil.copyfile(produced, dest)
    print(dest)

if __name__ == "__main__":
    main()
