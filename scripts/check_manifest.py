"""Check the distributed reader tree without requiring Git or a clean checkout."""
import argparse
import hashlib
from pathlib import Path, PurePosixPath
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = "MANIFEST_SHA256.txt"
IGNORED = {".git", "__pycache__", ".pytest_cache", "build", "venv", ".venv", "replay_work"}


def verify(root=ROOT):
    root = Path(root).resolve()
    entries = {}
    for row in (root/MANIFEST).read_text(encoding="utf8").splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", row)
        if match is None:
            raise ValueError("malformed manifest row")
        digest, name = match.groups()
        rel = PurePosixPath(name)
        if (rel.is_absolute() or ".." in rel.parts or "\\" in name or ":" in name
                or str(rel) != name or name == MANIFEST or name in entries
                or set(rel.parts) & IGNORED):
            raise ValueError("unsafe, duplicate or excluded manifest path: " + name)
        entries[name] = digest
    if not entries:
        raise ValueError("empty manifest")
    found = set()
    for p in root.rglob("*"):
        rel = p.relative_to(root)
        if set(rel.parts) & IGNORED:
            continue
        if p.is_symlink():
            raise ValueError("symlink in reader tree: " + rel.as_posix())
        if p.is_file():
            found.add(rel.as_posix())
    expected = set(entries) | {MANIFEST}
    if found != expected:
        raise ValueError("file-set mismatch: missing=" + repr(sorted(expected-found))
                         + ", unexpected=" + repr(sorted(found-expected)))
    for name, digest in entries.items():
        p = root/name
        if not p.resolve().is_relative_to(root):
            raise ValueError("path leaves reader tree")
        if hashlib.sha256(p.read_bytes()).hexdigest() != digest:
            raise ValueError("hash mismatch: " + name)
    return {"status": "MANIFEST_PASS", "files": len(found)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    try:
        import json
        print(json.dumps(verify(args.root), sort_keys=True))
    except (OSError, ValueError) as exc:
        print("MANIFEST_FAIL: " + str(exc), file=sys.stderr)
        return 1
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
