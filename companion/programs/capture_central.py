"""Run an unchanged aggregate producer while preserving its central arrays."""
import argparse
import hashlib
import importlib.util
import json
from math import isqrt
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("script", type=Path)
    parser.add_argument("source_sha256")
    parser.add_argument("batch", type=int, choices=(29, 30, 31, 32, 33, 35))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise ValueError("fresh central capture output required")
    if hashlib.sha256(args.script.read_bytes()).hexdigest() != args.source_sha256:
        raise ValueError("producer adapter byte identity mismatch")
    spec = importlib.util.spec_from_file_location("captured_aggregate", args.script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    original = module.central
    rows = {}

    def capture(*values, **keywords):
        result = original(*values, **keywords)
        M = values[0]
        p, q = ((2, 3 if args.batch == 29 else 5) if args.batch in (29, 30)
                else (values[1], values[2]))
        profile = [M, p, q]
        strata, stats = result
        row = {"M": M, "p": p, "r": q, "central_k": isqrt(M),
               "central_strata": [[g, list(cdf)] for g, cdf in sorted(strata.items())]}
        key = tuple(profile)
        if key in rows and rows[key] != row:
            raise ValueError("inconsistent repeated producer central call")
        rows[key] = row
        return result  # Exactly the original returned objects, without mutation.

    module.central = capture
    module.main()
    document = {"status": "PRODUCER_CENTRAL_ARRAY_CAPTURE_COMPLETE",
                "batch": args.batch, "producer_adapter_sha256": args.source_sha256,
                "rows": [rows[key] for key in sorted(rows)]}
    with args.output.open("x", encoding="utf8", newline="\n") as stream:
        json.dump(document, stream, sort_keys=True, separators=(",", ":"))
        stream.write("\n")


if __name__ == "__main__":
    main()
