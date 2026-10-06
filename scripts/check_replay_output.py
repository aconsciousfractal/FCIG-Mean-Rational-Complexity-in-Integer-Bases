"""Compare saved JSON outputs with the observed canonical replay identities."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def strict_json(text):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('duplicate JSON key: ' + key)
            result[key] = value
        return result
    def constant(value):
        raise ValueError('nonfinite JSON constant: ' + value)
    return json.loads(text, object_pairs_hook=pairs, parse_constant=constant)


def canonical_digest(path):
    value = strict_json(Path(path).read_text(encoding="utf8"))
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False).encode("utf8")
    return hashlib.sha256(encoded).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name", choices=("smoke", "retained", "semiprime", "oracle"))
    parser.add_argument("outputs", type=Path, nargs="+")
    args = parser.parse_args()
    expected = strict_json((ROOT / "companion/expected_replay.json").read_text(encoding="utf8"))
    digest = expected["outputs"][args.name]["canonical_json_sha256"]
    for path in args.outputs:
        actual = canonical_digest(path)
        if actual != digest:
            raise SystemExit("REPLAY_OUTPUT_FAIL: " + str(path))
    print(json.dumps({"status": "REPLAY_OUTPUT_PASS", "name": args.name,
                      "outputs": len(args.outputs)}, sort_keys=True))


if __name__ == "__main__":
    main()
