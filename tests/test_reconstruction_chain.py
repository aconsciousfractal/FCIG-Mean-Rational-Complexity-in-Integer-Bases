"""Lightweight control fixtures; no mathematical producer is executed."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import types
import unittest
from unittest import mock

PROGRAM = Path(__file__).resolve().parents[1] / "companion/programs/reconstruction_chain.py"
spec = importlib.util.spec_from_file_location("reconstruction_control_fixture", PROGRAM)
control = importlib.util.module_from_spec(spec)
spec.loader.exec_module(control)  # Standard-library-only control module.


class ChainFixtures(unittest.TestCase):
    def test_type_and_nonfinite_controls(self):
        self.assertFalse(control.same_types({"value": 1}, {"value": True}))
        self.assertFalse(control.same_types([1], [1.0]))
        for raw in (b'{"a":1,"a":2}', b'{"a":NaN}', b'{"a":1e999}'):
            with self.assertRaises(ValueError):
                control.strict(raw)

    def test_fresh_paths_and_newline_transport(self):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory)
            for name in ("../escape", "a/b", "", "a b", "C:drive"):
                with self.assertRaises(ValueError):
                    control.fresh_workspace(package, name)
            target = control.fresh_workspace(package, "good")
            target.mkdir(parents=True)
            with self.assertRaises(FileExistsError):
                control.fresh_workspace(package, "good")
        raw = b'{\r\n  "text":"literal\\n", "value":1\r\n}\r\n'
        actual = control.newline_transport(raw, b'{}\n')
        self.assertEqual(control.strict(actual), {"text": "literal\n", "value": 1})
        self.assertNotIn(b"\r\n", actual)

    def fixture(self, root):
        public = root / "public"
        output = root / "output"
        work = root / "work"
        public.mkdir()
        output.mkdir()
        work.mkdir()
        producer = "import pathlib,json\npathlib.Path('datum.json').write_bytes(b'{\\\"value\\\":7}\\n')\npathlib.Path('generated.marker').write_text('fresh')\nprint(json.dumps({'value':7},separators=(',',':')))\n"
        consumer = "import pathlib,json\nif pathlib.Path('generated.marker').read_text()!='fresh': raise ValueError('not generated')\nif json.loads(pathlib.Path('datum.json').read_bytes())!={'value':7}: raise ValueError('stale datum')\nprint(json.dumps({'read_generated':True},separators=(',',':')))\n"
        texts = {"producer.py": producer.encode(), "consumer.py": consumer.encode()}
        expected = {"datum.json": b'{"value":7}\n', "consumer.json": b'{"read_generated":true}\n'}
        jobs = []
        for name, target, dependencies in (("producer", "datum.json", []), ("consumer", "consumer.json", ["producer"])):
            source = name + ".py"
            pin = hashlib.sha256(texts[source]).hexdigest()
            jobs.append({"id": name, "batch": 0, "source_locator": source,
                "public_source_sha256": pin, "adapter_sha256": pin,
                "requires": dependencies,
                "inputs": [] if name == "producer" else [{"locator": "datum.json", "sha256": hashlib.sha256(expected["datum.json"]).hexdigest()}],
                "direct_argv": ["python", "-B", "${workspace}/_portable_" + source],
                "output": {"consumed_locator": target, "emission": "stdout_capture",
                           "retained_sha256": hashlib.sha256(expected[target]).hexdigest()}})
        plan = root / "plan.json"
        plan.write_text(json.dumps({"schema": "FULL_RECONSTRUCTION_PLAN_V1", "jobs": jobs}), encoding="utf8")
        contract = {"public_root": str(public), "output_root": str(output), "work_root": str(work),
            "plan": str(plan), "plan_sha256": control.sha(plan), "controller_sha256": control.sha(PROGRAM),
            "source_pins": {}, "jobs": [j["id"] for j in jobs], "modes": ["normal", "optimized"],
            "runtime_contract_sha256": "0" * 64, "per_job_timeout_seconds": 10,
            "max_output_bytes": 10000, "max_stderr_bytes": 10000}
        archive = types.SimpleNamespace(manifest=lambda: {"entries": [{"locator": name} for name in expected]},
                                        read_object=lambda record: expected[record["locator"]])
        def prepare(c, workspace, mode):
            if workspace.exists():
                return
            workspace.mkdir()
            for name, content in texts.items():
                (workspace / name).write_bytes(content)
                (workspace / ("_portable_" + name)).write_bytes(content)
            # Deliberately stale fixture datum: only producer can replace it.
            (workspace / "datum.json").write_bytes(b'{"value":11}\n')
            (workspace / "consumer.json").write_bytes(expected["consumer.json"])
        return contract, archive, prepare, output, work

    def test_generator_output_feeds_consumer_and_checkpoint_rejects_mutation(self):
        with tempfile.TemporaryDirectory() as directory:
            contract, archive, prepare, output, work = self.fixture(Path(directory))
            workspace, rows = control.run_jobs(contract, "normal", archive=archive, prepare=prepare)
            self.assertEqual([row["status"] for row in rows], ["PASS", "PASS"])
            self.assertEqual(control.strict((workspace / "datum.json").read_bytes()), {"value": 7})
            (workspace / "datum.json").write_bytes(b'{"value":8}\n')
            with self.assertRaisesRegex(ValueError, "checkpoint produced bytes"):
                control.run_jobs(contract, "normal", archive=archive, prepare=prepare)

    def test_bad_ids_prerequisite_mode_and_changed_plan_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            contract, archive, prepare, output, work = self.fixture(Path(directory))
            for ids in (["missing"], ["producer", "producer"], ["consumer"]):
                with self.assertRaises(ValueError):
                    control.run_jobs(contract, "normal", ids=ids, archive=archive, prepare=prepare)
            with self.assertRaises(ValueError):
                control.run_jobs(contract, "other", archive=archive, prepare=prepare)
            Path(contract["plan"]).write_text("{}", encoding="utf8")
            with self.assertRaisesRegex(ValueError, "byte freeze"):
                control.run_jobs(contract, "normal", archive=archive, prepare=prepare)

    def test_source_pin_and_receipt_mode_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            contract, archive, prepare, output, work = self.fixture(Path(directory))
            source = Path(contract["public_root"]) / "pinned.txt"
            source.write_bytes(b"original")
            contract["source_pins"] = {"pinned.txt": control.sha(source)}
            source.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "public source byte freeze"):
                control.run_jobs(contract, "normal", archive=archive, prepare=prepare)
            contract["source_pins"] = {}
            control.run_jobs(contract, "normal", archive=archive, prepare=prepare)
            receipt = output / "chain_normal_receipt.json"
            data = json.loads(receipt.read_bytes())
            data["mode"] = "optimized"
            receipt.write_text(json.dumps(data), encoding="utf8")
            with self.assertRaisesRegex(ValueError, "receipt contract/mode"):
                control.run_jobs(contract, "normal", archive=archive, prepare=prepare)

    def test_complete_status_requires_all_jobs(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "all 37"):
                control.complete_stages({"output_root": directory}, "normal", Path(directory), [], "full")

    def test_inventory_hook_derives_then_feeds_resolver(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            public, output, work = (root / name for name in ("public", "output", "work"))
            for folder in (public, output, work):
                folder.mkdir()
            programs = public / "companion/programs"
            programs.mkdir(parents=True)
            helper = programs / "derive_inventory_certificate.py"
            helper.write_text("# mocked metadata derivation\n", encoding="utf8")
            inventory = b'{"rows":[1]}\n'
            inventory_pin = hashlib.sha256(inventory).hexdigest()
            frontier = {"status": "R36_TWO_METHOD_DOMAIN_PASS", "source": "generated",
                        "generator_sha256": inventory_pin, "factor_sha256": inventory_pin}
            frontier_raw = (json.dumps(frontier, sort_keys=True, indent=2) + "\n").encode()
            base = "results/batch_36/"
            targets = {base + "inventory_generate_normal1.json": inventory,
                       base + "inventory_factor_normal1.json": inventory,
                       base + "FRONTIER_CERTIFICATE.json": frontier_raw,
                       base + "resolve_work_normal1.json": b'{"fresh_frontier":true}\n'}
            source_texts = {}
            jobs = []
            for name, filename in (("inventory_generate", "inventory_generate_normal1.json"),
                                   ("inventory_factor", "inventory_factor_normal1.json"),
                                   ("resolve_work36", "resolve_work_normal1.json")):
                if name == "resolve_work36":
                    text = "import pathlib,json\np=json.loads(pathlib.Path('results/batch_36/FRONTIER_CERTIFICATE.json').read_bytes())\nif p.get('source')!='generated': raise ValueError('stale frontier')\nprint(json.dumps({'fresh_frontier':True}))\n"
                    prerequisites = ["inventory_generate", "inventory_factor"]
                    inputs = [{"locator": base + "FRONTIER_CERTIFICATE.json", "sha256": hashlib.sha256(frontier_raw).hexdigest()}]
                else:
                    text = "import sys\nsys.stdout.write('{\\\"rows\\\":[1]}\\n')\n"
                    prerequisites, inputs = [], []
                locator = base + name + ".py"
                source_texts[locator] = text.encode()
                pin = hashlib.sha256(text.encode()).hexdigest()
                jobs.append({"id": name, "batch": 36, "source_locator": locator,
                    "public_source_sha256": pin, "adapter_sha256": pin, "requires": prerequisites,
                    "inputs": inputs, "direct_argv": ["python", "-B", "${workspace}/" + base + "_portable_" + name + ".py"],
                    "output": {"consumed_locator": base + filename, "emission": "stdout_capture",
                               "retained_sha256": hashlib.sha256(targets[base + filename]).hexdigest()}})
            plan = root / "plan.json"
            plan.write_text(json.dumps({"schema": "FULL_RECONSTRUCTION_PLAN_V1", "jobs": jobs}), encoding="utf8")
            contract = {"public_root": str(public), "output_root": str(output), "work_root": str(work),
                "plan": str(plan), "plan_sha256": control.sha(plan), "controller_sha256": control.sha(PROGRAM),
                "source_pins": {"companion/programs/derive_inventory_certificate.py": control.sha(helper)},
                "jobs": [j["id"] for j in jobs], "modes": ["normal"], "runtime_contract_sha256": "0" * 64,
                "per_job_timeout_seconds": 10, "max_output_bytes": 10000, "max_stderr_bytes": 10000}
            archive = types.SimpleNamespace(manifest=lambda: {"entries": [{"locator": n} for n in targets]},
                                            read_object=lambda record: targets[record["locator"]])
            def prepare(c, workspace, mode):
                (workspace / base).mkdir(parents=True)
                for locator, raw in source_texts.items():
                    target = workspace / locator
                    target.write_bytes(raw)
                    target.with_name("_portable_" + target.name).write_bytes(raw)
                for locator, raw in targets.items():
                    (workspace / locator).write_bytes(raw)
                (workspace / base / "FRONTIER_CERTIFICATE.json").write_bytes(b'{"source":"stale"}\n')
            real_execute = control.execute
            calls = []
            def execute(argv, workspace, stdout, stderr, c):
                if str(helper) not in argv:
                    return real_execute(argv, workspace, stdout, stderr, c)
                calls.append(argv)
                generated = Path(argv[-1])
                generated.write_bytes(frontier_raw)
                stdout.write_text(json.dumps({"certificate_sha256": control.sha(generated)}), encoding="utf8")
                stderr.write_bytes(b"")
                return {"exit_code": 0, "timeout_or_output_limit": False}
            with mock.patch.object(control, "execute", side_effect=execute):
                workspace, rows = control.run_jobs(contract, "normal", archive=archive, prepare=prepare)
            self.assertEqual(len(calls), 1)
            self.assertEqual([r["status"] for r in rows], ["PASS"] * 3)
            self.assertEqual((workspace / base / "FRONTIER_CERTIFICATE.json").read_bytes(), frontier_raw)
            self.assertEqual(control.strict((workspace / base / "resolve_work_normal1.json").read_bytes()), {"fresh_frontier": True})

    def stage_fixture(self, root, fault=None):
        workspace, output = root / "workspace", root / "output"
        workspace.mkdir()
        output.mkdir()
        counts = {29:143, 30:97, 31:132, 32:3, 33:35, 35:99, 39:1087, 40:1422, 41:1901}
        jobs, rows = [], []
        for name in sorted(control.REQUIRED_JOB_IDS):
            batch = int(name[len("aggregate"):]) if name.startswith("aggregate") else 0
            locator = name + ".json"
            target = workspace / locator
            target.write_bytes(b'{"fixture":true}\n')
            row = {"id": name, "status": "PASS", "attempt": 1,
                   "consumed_output_sha256": control.sha(target)}
            if batch in control.EARLY_BATCHES:
                capture = output / (name + "_normal_central.json")
                capture.write_bytes(b'{}\n')
                row["capture_sha256"] = control.sha(capture)
            rows.append(row)
            jobs.append({"id": name, "source_locator": name + ".py", "inputs": [], "requires": [],
                         "scope_count": counts.get(batch, 1), "output": {"consumed_locator": locator}})
        for name, filename in (("base6", "results/batch_26/FINITE_BRIDGE.json"),
                               ("sixth", "results/batch_27/SIXTH_LENGTH_BRIDGE.json")):
            target = workspace / filename
            target.parent.mkdir(parents=True)
            target.write_bytes(b'{"fixture":true}\n')
            next(row for row in rows if row["id"] == name)["consumed_output_sha256"] = control.sha(target)
        plan = root / "plan.json"
        plan.write_text(json.dumps({"schema": "FULL_RECONSTRUCTION_PLAN_V1", "jobs": jobs}), encoding="utf8")
        contract = {"output_root": str(output), "plan": str(plan)}
        calls = []
        def external(c, mode, work, out, name, program, args):
            calls.append(name)
            if program == "reconstruct_semiprime.py":
                path = Path(args[-1])
                path.write_bytes(b'{"fixture":"semiprime"}\n')
                summary = {"status": "PORTABLE_E032_RECONSTRUCTION_PASS", "scope": "full",
                    "reconstructed_rows": 4596, "output_sha256": control.sha(path), "output_bytes": path.stat().st_size}
                if fault == "semiprime_pin":
                    summary["output_sha256"] = "0" * 64
            elif program == "verify_full_semiprime.py":
                summary = {"status": "COMPLETE_SEMIPRIME_HISTOGRAM_AND_DISTINCT_LOG_PASS",
                    "histograms_checked": 4596, "input_sha256": args[1],
                    "upper_envelope_entries": 41585055, "every_outward_sum_negative": True}
            elif program == "verify_all_central.py":
                batch = int(args[args.index("--batch") + 1])
                count = counts[batch]
                report = {"status": "COMPLETE_REQUESTED_CHUNK_PASS", "batch": batch,
                    "complete_batch_verified": True, "complete_batch_rows": count, "start": 0, "stop": count,
                    "coverage_indices": list(range(count)), "rows": [{"domain_index": i} for i in range(count)],
                    "pins": {"consumed_certificate_sha256": args[args.index("--certificate-sha256") + 1]}}
                if fault == "central_flag":
                    report["complete_batch_verified"] = False
                if fault == "central_missing":
                    report["rows"].pop()
                Path(args[args.index("--output") + 1]).write_text(json.dumps(report), encoding="utf8")
                summary = {"status": report["status"]}
            else:
                kind = args[args.index("--kind") + 1]
                count = {"semiprime":4596, "base6":8, "sixth":3}[kind]
                report = {"status": "COMPLETE_POINTWISE_CENTRAL_CHUNK_PASS", "kind": kind,
                    "complete_domain_verified": True, "full_domain_count": count, "start":0, "stop":count,
                    "rows": [{"domain_index": i} for i in range(count)],
                    "input_sha256": args[args.index("--input-sha256") + 1]}
                if fault == "pointwise_pin":
                    report["input_sha256"] = "0" * 64
                Path(args[args.index("--output") + 1]).write_text(json.dumps(report), encoding="utf8")
                summary = {"status": report["status"]}
            return {"name": name, "status": "PASS", "summary": summary}
        return contract, workspace, rows, external, calls

    def test_full_and_partial_stage_claims(self):
        for stage, complete in (("full", True), ("reconstruction", False), ("central", False)):
            with self.subTest(stage=stage), tempfile.TemporaryDirectory() as directory:
                contract, workspace, rows, external, calls = self.stage_fixture(Path(directory))
                with mock.patch.object(control, "external", side_effect=external):
                    result = control.complete_stages(contract, "normal", workspace, rows, stage)
                self.assertIs(result["complete"], complete)
                self.assertEqual("semiprime_distinct_logs" in calls, stage == "full")
                self.assertEqual(sum(name.startswith("central_") for name in calls), 0 if stage == "reconstruction" else 9)
                self.assertEqual(sum(name.startswith("pointwise_") for name in calls), 0 if stage == "reconstruction" else 3)

    def test_incomplete_or_unbound_distinct_reports_never_complete(self):
        for fault in ("semiprime_pin", "central_flag", "central_missing", "pointwise_pin"):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as directory:
                contract, workspace, rows, external, calls = self.stage_fixture(Path(directory), fault)
                with mock.patch.object(control, "external", side_effect=external):
                    with self.assertRaises(ValueError):
                        control.complete_stages(contract, "normal", workspace, rows, "full")


if __name__ == "__main__":
    unittest.main()
