"""Portable reconstruction control derived from the qualified v3 controller.

Producer arithmetic stays in the hash-pinned public adapters. Only generated
work contracts and receipts contain the reader's current filesystem paths.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import time

V3_CONTROL_SHA256 = "5a66a18e80ea1ab81009551d9fc627c0935300ea9bcad52ea5459bbeefc111b5"
EARLY_BATCHES = (29, 30, 31, 32, 33, 35)
ALL_BATCHES = EARLY_BATCHES + (39, 40, 41)
REQUIRED_JOB_IDS = frozenset((
    "prime3", "prime13", "mixed_fixed", "counterexample", "base6", "sixth", "semiprime_all",
    "aggregate29", "aggregate30", "aggregate31", "aggregate32", "aggregate33",
    "squares_analytic", "squares_check", "aggregate35", "inventory_generate", "inventory_factor",
    "resolve_work36", "n1_analytic", "n1_check", "remaining_work37", "resolve_cost38",
    "aggregate39", "check39", "physical39", "aggregate40", "check40", "physical40",
    "aggregate41", "check41", "physical41", "retained_check29", "retained_check30",
    "retained_check31", "retained_check32", "retained_check33", "retained_check35"))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            digest.update(chunk)
    return digest.hexdigest()


def strict(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("duplicate JSON key: " + key)
            result[key] = value
        return result
    def constant(value):
        raise ValueError("nonfinite JSON constant: " + value)
    def finite_float(value):
        number = float(value)
        if not math.isfinite(number):
            raise ValueError("nonfinite JSON number")
        return number
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=constant,
                      parse_float=finite_float)


def same_types(left, right):
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return set(left) == set(right) and all(same_types(left[k], right[k]) for k in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(same_types(a, b) for a, b in zip(left, right))
    return left == right


def relative(name):
    if type(name) is not str:
        raise ValueError("relative locator string required")
    path = PurePosixPath(name)
    if (not name or path.is_absolute() or ".." in path.parts or "\\" in name
            or ":" in name or str(path) != name or name == "."):
        raise ValueError("unsafe relative locator: " + name)
    return name


def fresh_workspace(package, name):
    if type(name) is not str or re.fullmatch(r"[A-Za-z0-9_-]+", name) is None:
        raise ValueError("work name must contain only letters, digits, underscore or hyphen")
    package = Path(package).resolve()
    parent = package / "replay_work"
    target = parent / name
    if parent.is_symlink() or target.exists() or target.is_symlink():
        raise FileExistsError("a fresh nonsymlink work directory is required")
    if not target.resolve().is_relative_to(package):
        raise ValueError("work directory escapes package")
    return target


def newline_transport(raw, baseline):
    lf = raw.replace(b"\r\n", b"\n")
    transported = lf.replace(b"\n", b"\r\n") if b"\r\n" in baseline else lf
    if not same_types(strict(raw), strict(transported)):
        raise ValueError("newline transport changed JSON meaning")
    return transported


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n", encoding="utf8", newline="\n")


def load_plan(path):
    plan = strict(Path(path).read_bytes())
    if plan.get("schema") != "FULL_RECONSTRUCTION_PLAN_V1" or type(plan.get("jobs")) is not list:
        raise ValueError("reconstruction plan schema")
    ids = set()
    for job in plan["jobs"]:
        if type(job.get("id")) is not str or job["id"] in ids:
            raise ValueError("duplicate or invalid job ID")
        ids.add(job["id"])
        relative(job["source_locator"])
        relative(job["output"]["consumed_locator"])
        for record in job["inputs"]:
            relative(record["locator"])
        if type(job["requires"]) is not list:
            raise ValueError("job prerequisites must be a list")
    if any(not set(job["requires"]) <= ids for job in plan["jobs"]):
        raise ValueError("unknown prerequisite ID")
    return plan


def archive_modules(public_root):
    sys.path.insert(0, str(Path(public_root) / "companion" / "programs"))
    import archive
    import replay
    return archive, replay


def materialize(contract, workspace, mode, archive, replay):
    origin_path = workspace / "SOURCE_ORIGIN.json"
    contract_pin = contract["runtime_contract_sha256"]
    if workspace.exists():
        origin = strict(origin_path.read_bytes())
        if origin.get("contract_sha256") != contract_pin or origin.get("mode") != mode:
            raise ValueError("work origin contract/mode mismatch")
        return
    objects, adapters = replay.preflight(workspace, archive.manifest()["entries"],
        strict((Path(contract["public_root"]) / "companion/portable_source_map.json").read_bytes()))
    workspace.mkdir(parents=True, exist_ok=False)
    for target, record in objects:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(archive.read_object(record))
    for target, payload in adapters:
        target.write_bytes(payload)
    write_json(origin_path, {"contract_sha256": contract_pin, "mode": mode})


def execute(argv, workspace, output, error, contract):
    if output.exists() or error.exists():
        raise FileExistsError("fresh attempt output required")
    started = time.monotonic()
    exceeded = False
    with output.open("xb") as stdout, error.open("xb") as stderr:
        process = subprocess.Popen(argv, cwd=workspace, stdout=stdout, stderr=stderr)
        while process.poll() is None:
            if (time.monotonic() - started > contract["per_job_timeout_seconds"]
                    or output.stat().st_size > contract["max_output_bytes"]
                    or error.stat().st_size > contract["max_stderr_bytes"]):
                process.kill()
                exceeded = True
                break
            time.sleep(.25)
        code = process.wait()
    exceeded = exceeded or output.stat().st_size > contract["max_output_bytes"] or error.stat().st_size > contract["max_stderr_bytes"]
    return {"exit_code": code, "seconds": time.monotonic() - started,
            "timeout_or_output_limit": exceeded, "stdout_sha256": sha(output),
            "stdout_bytes": output.stat().st_size, "stderr_sha256": sha(error)}


def pinned_program(contract, name):
    locator = "companion/programs/" + name
    path = Path(contract["public_root"]) / locator
    if sha(path) != contract["source_pins"][locator]:
        raise ValueError("program byte pin changed: " + name)
    return path


def baseline_for(archive, locator):
    records = [r for r in archive.manifest()["entries"] if r["locator"] == locator]
    if len(records) != 1:
        raise ValueError("ambiguous output baseline: " + locator)
    return archive.read_object(records[0])


def derive_frontier(contract, workspace, output_root, mode, archive):
    destination = workspace / "results/batch_36/FRONTIER_CERTIFICATE.json"
    output = output_root / "inventory_certificate_generated.json"
    note = output_root / "inventory_certificate_receipt.json"
    generator = workspace / "results/batch_36/inventory_generate_normal1.json"
    factor = workspace / "results/batch_36/inventory_factor_normal1.json"
    program = pinned_program(contract, "derive_inventory_certificate.py")
    binding = {"generator_sha256": sha(generator), "factor_sha256": sha(factor),
               "program_sha256": sha(program), "mode": mode}
    if note.exists():
        receipt = strict(note.read_bytes())
        if any(receipt.get(k) != v for k, v in binding.items()) or receipt.get("status") != "PASS":
            raise ValueError("inventory derivation checkpoint binding changed")
        if sha(output) != receipt["generated_sha256"] or sha(destination) != receipt["consumed_sha256"]:
            raise ValueError("inventory derivation checkpoint bytes changed")
        return receipt
    argv = [sys.executable, "-B"] + (["-O"] if mode == "optimized" else []) + [
        str(program), str(generator), binding["generator_sha256"], str(factor),
        binding["factor_sha256"], str(output)]
    stdout = output_root / "inventory_certificate_summary.json"
    stderr = output_root / "inventory_certificate.stderr"
    receipt = binding | execute(argv, workspace, stdout, stderr, contract) | {"status": "FAILED"}
    try:
        if receipt["exit_code"] or receipt["timeout_or_output_limit"]:
            raise ValueError("inventory certificate derivation failed")
        summary = strict(stdout.read_bytes())
        if summary.get("certificate_sha256") != sha(output):
            raise ValueError("inventory certificate summary binding")
        baseline = baseline_for(archive, "results/batch_36/FRONTIER_CERTIFICATE.json")
        if not same_types(strict(output.read_bytes()), strict(baseline)):
            raise ValueError("derived inventory certificate differs from complete baseline")
        destination.write_bytes(newline_transport(output.read_bytes(), baseline))
        receipt.update(status="PASS", generated_sha256=sha(output), consumed_sha256=sha(destination))
    finally:
        write_json(note, receipt)
    return receipt


def run_jobs(contract, mode, ids=None, archive=None, prepare=None):
    if mode not in ("normal", "optimized") or mode not in contract["modes"]:
        raise ValueError("invalid reconstruction mode")
    public_root = Path(contract["public_root"])
    output_root = Path(contract["output_root"])
    workspace = Path(contract["work_root"]) / ("chain_" + mode)
    if sha(Path(__file__)) != contract["controller_sha256"] or sha(contract["plan"]) != contract["plan_sha256"]:
        raise ValueError("controller/plan byte freeze changed")
    for name, pin in contract["source_pins"].items():
        if sha(public_root / relative(name)) != pin:
            raise ValueError("public source byte freeze changed: " + name)
    jobs = load_plan(contract["plan"])["jobs"]
    by_id = {job["id"]: job for job in jobs}
    ids = list(contract["jobs"] if ids is None else ids)
    if len(set(ids)) != len(ids) or not set(ids) <= set(contract["jobs"]) <= set(by_id):
        raise ValueError("unknown/duplicate requested job")
    if archive is None:
        archive, replay = archive_modules(public_root)
    else:
        replay = None
    if prepare is None:
        materialize(contract, workspace, mode, archive, replay)
    else:
        prepare(contract, workspace, mode)
    receipt_path = output_root / ("chain_" + mode + "_receipt.json")
    rows = []
    if receipt_path.exists():
        prior = strict(receipt_path.read_bytes())
        if prior.get("mode") != mode or prior.get("contract_sha256") != contract["runtime_contract_sha256"]:
            raise ValueError("receipt contract/mode mismatch")
        rows = prior["jobs"]
    def target_for(job):
        locator = job["output"]["consumed_locator"]
        if job["id"] == "counterexample" and mode == "optimized":
            locator = locator.replace("_normal.json", "_O.json")
        return workspace / locator
    done = set()
    for row in rows:
        if row.get("status") != "PASS":
            continue
        if row.get("mode") != mode or row.get("id") not in by_id or row.get("controller_sha256") != contract["controller_sha256"]:
            raise ValueError("checkpoint row identity mismatch")
        job = by_id[row["id"]]
        suffix = "" if row.get("attempt", 1) == 1 else "_a" + str(row["attempt"])
        if sha(target_for(job)) != row["consumed_output_sha256"] or sha(output_root / (job["id"] + "_" + mode + suffix + ".json")) != row["stdout_sha256"]:
            raise ValueError("checkpoint produced bytes changed")
        if row.get("capture_sha256") and sha(output_root / (job["id"] + "_" + mode + suffix + "_central.json")) != row["capture_sha256"]:
            raise ValueError("checkpoint capture bytes changed")
        done.add(job["id"])
    for job in jobs:
        if job["id"] not in ids or job["id"] in done:
            continue
        if not set(job["requires"]) <= done:
            raise ValueError("uncompleted prerequisite: " + job["id"])
        if job["id"] == "resolve_work36":
            derive_frontier(contract, workspace, output_root, mode, archive)
        input_pins = {}
        for record in job["inputs"]:
            actual = sha(workspace / record["locator"])
            if actual != record["sha256"]:
                raise ValueError("actual consumed input byte drift: " + record["locator"])
            input_pins[record["locator"]] = actual
        original = workspace / job["source_locator"]
        adapter = original.with_name("_portable_" + original.name)
        if sha(original) != job["public_source_sha256"] or sha(adapter) != job["adapter_sha256"]:
            raise ValueError("materialized source byte drift")
        def argument(value):
            if value.startswith("${sha256:"):
                return sha(workspace / value[len("${sha256:"):-1])
            return value.replace("${workspace}", str(workspace))
        argv = [sys.executable, "-B"] + (["-O"] if mode == "optimized" else []) + [argument(v) for v in job["direct_argv"][2:]]
        attempt = 1 + sum(row["id"] == job["id"] for row in rows)
        suffix = "" if attempt == 1 else "_a" + str(attempt)
        capture = None
        if job["id"] in {"aggregate" + str(batch) for batch in EARLY_BATCHES}:
            capture = output_root / (job["id"] + "_" + mode + suffix + "_central.json")
            argv = [sys.executable, "-B"] + (["-O"] if mode == "optimized" else []) + [str(pinned_program(contract, "capture_central.py")), str(adapter), job["adapter_sha256"], str(job["batch"]), str(capture)]
        output = output_root / (job["id"] + "_" + mode + suffix + ".json")
        error = output_root / (job["id"] + "_" + mode + suffix + ".stderr")
        target = target_for(job)
        baseline = baseline_for(archive, target.relative_to(workspace).as_posix())
        print("START " + mode + " " + job["id"], flush=True)
        row = {"id": job["id"], "attempt": attempt, "mode": mode,
               "controller_sha256": sha(Path(__file__)), "argv": argv, "input_sha256": input_pins,
               "status": "FAILED"}
        try:
            row.update(execute(argv, workspace, output, error, contract))
            row["capture_sha256"] = sha(capture) if capture and capture.exists() else None
            if row["exit_code"] or row["timeout_or_output_limit"]:
                raise ValueError("producer job failed: " + job["id"])
            summary = strict(output.read_bytes())
            rebuilt = strict(target.read_bytes()) if job["id"] == "counterexample" else summary
            if job["id"] == "counterexample" and (summary.get("sha256") != sha(target) or summary.get("output") != target.name):
                raise ValueError("summary/full-file binding")
            if not same_types(rebuilt, strict(baseline)):
                raise ValueError("retained/rebuilt complete value or type drift: " + job["id"])
            if job["output"]["emission"] == "stdout_capture":
                target.write_bytes(newline_transport(output.read_bytes(), baseline))
                row["stdout_transport"] = "newline whitespace only; raw stdout preserved"
            elif not same_types(strict(target.read_bytes()), rebuilt):
                raise ValueError("file/stdout disagreement")
            row.update(consumed_output_sha256=sha(target), retained_math_equal=True,
                       raw_retained_equal=sha(target) == job["output"]["retained_sha256"], status="PASS")
            done.add(job["id"])
        except Exception as failure:
            row["error"] = {"type": type(failure).__name__, "message": str(failure)}
            raise
        finally:
            rows.append(row)
            write_json(receipt_path, {"contract_sha256": contract["runtime_contract_sha256"], "mode": mode, "jobs": rows})
        print("PASS " + mode + " " + job["id"], flush=True)
    if not set(ids) <= done:
        raise ValueError("requested completion set incomplete")
    return workspace, rows


def external(contract, mode, workspace, output_root, name, program, args):
    argv = [sys.executable, "-B"] + (["-O"] if mode == "optimized" else []) + [str(pinned_program(contract, program)), *map(str, args)]
    stdout = output_root / (name + "_summary.json")
    stderr = output_root / (name + ".stderr")
    row = {"name": name, "argv": argv, "status": "FAILED"}
    try:
        row.update(execute(argv, workspace, stdout, stderr, contract))
        if row["exit_code"] or row["timeout_or_output_limit"]:
            raise ValueError("distinct stage failed: " + name)
        row["summary"] = strict(stdout.read_bytes())
        row["status"] = "PASS"
        return row
    finally:
        write_json(output_root / (name + "_receipt.json"), row)


def complete_stages(contract, mode, workspace, chain_rows, stage):
    if stage not in ("full", "reconstruction", "central"):
        raise ValueError("unknown stage")
    output_root = Path(contract["output_root"])
    by_id = {row["id"]: row for row in chain_rows if row.get("status") == "PASS"}
    if set(by_id) != REQUIRED_JOB_IDS:
        raise ValueError("all 37 reconstruction jobs are required before further stages")
    semiprime = output_root / "semiprime_full.json"
    reconstruct = external(contract, mode, workspace, output_root, "semiprime_full", "reconstruct_semiprime.py", ["--all", "--output", semiprime])
    summary = reconstruct["summary"]
    if (summary.get("status") != "PORTABLE_E032_RECONSTRUCTION_PASS" or summary.get("scope") != "full"
            or type(summary.get("reconstructed_rows")) is not int or summary["reconstructed_rows"] != 4596
            or summary.get("output_sha256") != sha(semiprime) or summary.get("output_bytes") != semiprime.stat().st_size):
        raise ValueError("complete semiprime reconstruction binding")
    stages = [reconstruct]
    if stage == "full":
        check = external(contract, mode, workspace, output_root, "semiprime_distinct_logs", "verify_full_semiprime.py", [semiprime, sha(semiprime)])
        if (check["summary"].get("status") != "COMPLETE_SEMIPRIME_HISTOGRAM_AND_DISTINCT_LOG_PASS"
                or check["summary"].get("histograms_checked") != 4596 or check["summary"].get("input_sha256") != sha(semiprime)
                or check["summary"].get("upper_envelope_entries") != 41585055 or check["summary"].get("every_outward_sum_negative") is not True):
            raise ValueError("complete distinct semiprime log gate")
        stages.append(check)
    if stage in ("full", "central"):
        for batch in ALL_BATCHES:
            producer = by_id["aggregate" + str(batch)]
            job = next(job for job in load_plan(contract["plan"])["jobs"] if job["id"] == producer["id"])
            locator = job["output"]["consumed_locator"]
            certificate = workspace / locator
            if sha(certificate) != producer["consumed_output_sha256"]:
                raise ValueError("fresh aggregate consumer binding changed")
            report = output_root / ("central_" + str(batch) + ".json")
            args = ["--batch", batch, "--certificate", certificate, "--certificate-sha256", sha(certificate), "--output", report]
            if batch in EARLY_BATCHES:
                suffix = "" if producer["attempt"] == 1 else "_a" + str(producer["attempt"])
                capture = output_root / (producer["id"] + "_" + mode + suffix + "_central.json")
                if sha(capture) != producer["capture_sha256"]:
                    raise ValueError("fresh central capture consumer binding changed")
                args += ["--producer-arrays", capture, "--producer-arrays-sha256", sha(capture)]
            checked = external(contract, mode, workspace, output_root, "central_" + str(batch), "verify_all_central.py", args)
            value = strict(report.read_bytes())
            if (value.get("status") != "COMPLETE_REQUESTED_CHUNK_PASS" or value.get("complete_batch_verified") is not True
                    or value.get("batch") != batch or value.get("pins", {}).get("consumed_certificate_sha256") != sha(certificate)
                    or value.get("complete_batch_rows") != job["scope_count"] or value.get("start") != 0 or value.get("stop") != job["scope_count"]
                    or value.get("coverage_indices") != list(range(job["scope_count"])) or type(value.get("rows")) is not list
                    or [row.get("domain_index") for row in value["rows"]] != list(range(job["scope_count"]))):
                raise ValueError("incomplete aggregate central domain")
            checked["report_sha256"] = sha(report)
            stages.append(checked)
        for kind, input_file, count in (("semiprime", semiprime, 4596), ("base6", workspace / "results/batch_26/FINITE_BRIDGE.json", 8), ("sixth", workspace / "results/batch_27/SIXTH_LENGTH_BRIDGE.json", 3)):
            expected_input = summary["output_sha256"] if kind == "semiprime" else by_id["base6" if kind == "base6" else "sixth"]["consumed_output_sha256"]
            if sha(input_file) != expected_input:
                raise ValueError("fresh pointwise consumer binding changed")
            report = output_root / ("pointwise_" + kind + ".json")
            checked = external(contract, mode, workspace, output_root, "pointwise_" + kind, "verify_pointwise_central.py", ["--kind", kind, "--input", input_file, "--input-sha256", sha(input_file), "--output", report])
            value = strict(report.read_bytes())
            if (value.get("status") != "COMPLETE_POINTWISE_CENTRAL_CHUNK_PASS" or value.get("complete_domain_verified") is not True
                    or value.get("full_domain_count") != count or value.get("kind") != kind or value.get("input_sha256") != expected_input
                    or value.get("start") != 0 or value.get("stop") != count or type(value.get("rows")) is not list
                    or [row.get("domain_index") for row in value["rows"]] != list(range(count))):
                raise ValueError("incomplete pointwise central domain")
            checked["report_sha256"] = sha(report)
            stages.append(checked)
    return {"status": "FULL_RECONSTRUCTION_AND_DISTINCT_CHECKS_PASS" if stage == "full" else "PARTIAL_RECONSTRUCTION_PASS",
            "complete": stage == "full", "stage": stage, "mode": mode,
            "chain_jobs": len(by_id), "stages": stages,
            "formal_proof_or_originality_certified": False}


def build_contract(public_root, work, mode):
    public_root = Path(public_root).resolve()
    plan = public_root / "companion/reconstruction_plan.json"
    parsed = load_plan(plan)
    if {job["id"] for job in parsed["jobs"]} != REQUIRED_JOB_IDS:
        raise ValueError("complete 37-job reconstruction plan required")
    names = ("archive.py", "replay.py", "capture_central.py", "derive_inventory_certificate.py",
             "reconstruct_semiprime.py", "sixth_length_math.py", "verify_full_semiprime.py",
             "verify_all_central.py", "verify_pointwise_central.py", "verify_retained.py",
             "reconstruction_chain.py")
    pins = {"companion/programs/" + name: sha(public_root / "companion/programs" / name) for name in names}
    pins.update({name: sha(public_root / name) for name in ("companion/provenance.json", "companion/portable_source_map.json")})
    return {"public_root": str(public_root), "work_root": str(work), "output_root": str(Path(work) / "output"),
            "plan": str(plan), "plan_sha256": sha(plan), "controller_sha256": sha(Path(__file__)),
            "source_pins": pins, "jobs": [job["id"] for job in parsed["jobs"]], "modes": [mode],
            "per_job_timeout_seconds": 600, "max_output_bytes": 128 * 1024 * 1024,
            "max_stderr_bytes": 1024 * 1024, "RAM_measured": False,
            "json_capture_transport": "newline whitespace only; raw stdout preserved"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", required=True, help="fresh simple package-local work name")
    parser.add_argument("--optimized", action="store_true")
    parser.add_argument("--stage", choices=("full", "reconstruction", "central"), default="full",
                        help="partial stages regenerate prerequisites; only full reports complete")
    args = parser.parse_args()
    public_root = Path(__file__).resolve().parents[2]
    work = fresh_workspace(public_root / "companion", args.work)
    mode = "optimized" if args.optimized else "normal"
    # Bind the initially distributed tree before constructing current-path pins.
    subprocess.run([sys.executable, "-B", str(public_root / "scripts/check_manifest.py")],
                   cwd=public_root, check=True)
    contract = build_contract(public_root, work, mode)
    work.mkdir(parents=True, exist_ok=False)
    Path(contract["output_root"]).mkdir()
    contract_path = work / "runtime_contract.json"
    write_json(contract_path, contract)
    contract["runtime_contract_sha256"] = sha(contract_path)
    final_path = Path(contract["output_root"]) / "completion.json"
    try:
        workspace, rows = run_jobs(contract, mode)
        result = complete_stages(contract, mode, workspace, rows, args.stage)
    except Exception as failure:
        write_json(final_path, {"status": "FAILED_INCOMPLETE_RECONSTRUCTION", "complete": False,
                              "mode": mode, "stage": args.stage,
                              "error": {"type": type(failure).__name__, "message": str(failure)}})
        raise
    write_json(final_path, result)
    print(json.dumps({key: result[key] for key in ("status", "complete", "stage", "mode", "chain_jobs")}, sort_keys=True))


if __name__ == "__main__":
    main()
