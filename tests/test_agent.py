"""Agent layer: status scanner, tools, number grounding, and the LangGraph loop.

The LLM is replaced by a scripted chat model, so these tests need no API key and make no network calls.
"""

import json

import pytest

pytest.importorskip("langgraph")

from langchain_core.language_models.chat_models import BaseChatModel  # noqa: E402
from langchain_core.messages import AIMessage, HumanMessage  # noqa: E402
from langchain_core.outputs import ChatGeneration, ChatResult  # noqa: E402
from pymatgen.core import Lattice, Structure  # noqa: E402

from elastemp.agent.graph import ask, build_graph  # noqa: E402
from elastemp.agent.grounding import ungrounded_numbers  # noqa: E402
from elastemp.agent.status import project_status  # noqa: E402
from elastemp.agent.tools import tool_inspect_run, tool_run_step  # noqa: E402
from test_convergence import _vasprun  # noqa: E402

STRAINS = [-0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03]


# --------------------------------------------------------------------------- fixtures


def _project(tmp_path, *, vasp=True):
    """A cubic project after make_elastic, with (optionally) finished static VASP runs."""
    ws = tmp_path / "project"
    ws.mkdir()
    Structure(Lattice.cubic(3.57), ["C"], [[0, 0, 0]]).to(filename=str(ws / "POSCAR"), fmt="poscar")
    (ws / "POTCAR").write_text("dummy")
    (ws / "symmetry").write_text("cubic")
    (ws / "strains.txt").write_text("".join(f"{e}\n" for e in STRAINS))
    for d in ["deformation-1", "deformation-2", "deformation-3", "deformation-bulk"]:
        for i in range(1, len(STRAINS) + 1):
            folder = ws / d / f"strain-{i}"
            folder.mkdir(parents=True)
            if vasp:
                (folder / "vasprun.xml").write_text(_vasprun(11))
                (folder / "OSZICAR").write_text("   1 F= -.71E+02 E0= -.71E+02  d E =0.0E+00\n")
    return ws


class ScriptedModel(BaseChatModel):
    """Plays back a fixed list of AIMessages and records what it was sent."""

    responses: list
    received: list = []

    @property
    def _llm_type(self):
        return "scripted"

    def bind_tools(self, tools, **kwargs):
        return self

    def _generate(self, messages, stop=None, run_manager=None, **kwargs):
        self.received.append(list(messages))
        return ChatResult(generations=[ChatGeneration(message=self.responses.pop(0))])


def _call(name, args=None, call_id="c1"):
    return AIMessage(content="", tool_calls=[{"name": name, "args": args or {}, "id": call_id}])


# --------------------------------------------------------------------------- status


def test_status_after_static_runs_says_run_get_elastic(tmp_path):
    report = project_status(_project(tmp_path))
    stages = {s["name"]: s["status"] for s in report["stages"]}
    assert report["symmetry"] == "cubic" and report["strains"] == 7
    assert stages["make_elastic"] == "done"
    assert stages["VASP static runs (deformation-*/strain-*)"] == "done"
    assert report["next_action"] == "Run the get_elastic step."


def test_status_lists_runs_still_to_do(tmp_path):
    ws = _project(tmp_path)
    for f in ("vasprun.xml", "OSZICAR"):
        (ws / "deformation-2" / "strain-5" / f).unlink()
    report = project_status(ws)
    static = report["stages"][1]
    assert static["status"] == "partial"
    assert static["not_run"] == ["deformation-2: strain-5"]
    assert report["next_action"].startswith("Run VASP")


def test_status_summarises_unrun_folders_completely(tmp_path):
    ws = _project(tmp_path, vasp=False)  # straight after make_elastic: nothing run yet
    static = project_status(ws)["stages"][1]
    assert static["detail"] == "0 of 28 converged, 0 failed, 28 not run"
    assert static["not_run"] == [
        "deformation-1: all 7 strains",
        "deformation-2: all 7 strains",
        "deformation-3: all 7 strains",
        "deformation-bulk: all 7 strains",
    ]


def test_status_blocks_when_a_deformation_has_too_few_points(tmp_path):
    ws = _project(tmp_path)  # 7 strains: any failure leaves fewer than 7 points
    (ws / "deformation-3" / "strain-2" / "vasprun.xml").write_text(_vasprun(60, nelm=60))
    static = project_status(ws)["stages"][1]
    assert static["status"] == "blocked"
    assert "deformation-3" in static["detail"]
    assert "NELM" in static["failed"][0]


# --------------------------------------------------------------------------- tools


def test_run_step_refuses_to_delete_finished_runs(tmp_path):
    ws = _project(tmp_path)
    result = json.loads(tool_run_step(ws, "make_elastic"))
    assert result["status"] == "refused"
    assert (ws / "deformation-1" / "strain-1" / "vasprun.xml").exists()


def test_run_step_rejects_unknown_steps(tmp_path):
    assert json.loads(tool_run_step(_project(tmp_path), "rm_rf"))["status"] == "invalid"


def test_inspect_run_stays_inside_project(tmp_path):
    ws = _project(tmp_path)
    assert "outside" in json.loads(tool_inspect_run(ws, "../.."))["error"]
    assert json.loads(tool_inspect_run(ws, "deformation-1/strain-1"))["converged"] is True


# --------------------------------------------------------------------------- grounding


def test_grounding_accepts_rounded_and_signless_values():
    refs = ['{"c11": 1054.086, "change": -4.54, "gap": 3.09}']
    assert ungrounded_numbers("C11 is 1054.1 GPa (about 1054 GPa); it drops 4.5%; gap 3.1%.", refs) == []


def test_grounding_ignores_identifiers_and_small_counts():
    assert ungrounded_numbers("deformation-3, strain-4 and C44: 3 runs failed in step 7.", []) == []


def test_grounding_catches_invented_values():
    assert ungrounded_numbers("Experimentally C11 = 1079 GPa.", ['{"c11": 1054.086}']) == ["1079"]


# --------------------------------------------------------------------------- graph


def test_agent_calls_tool_then_answers(tmp_path):
    ws = _project(tmp_path)
    model = ScriptedModel(
        responses=[_call("project_status"), AIMessage("All 28 static runs converged. Next: run get_elastic.")],
        received=[],
    )
    events = []
    answer = ask(build_graph(ws, model), "Where am I?", on_event=lambda k, d: events.append((k, d)))
    assert answer == "All 28 static runs converged. Next: run get_elastic."
    assert events == [("tool", "project_status()")]
    tool_output = model.received[1][-1].content
    assert '"next_action": "Run the get_elastic step."' in tool_output


def test_invented_number_triggers_one_rewrite(tmp_path):
    model = ScriptedModel(
        responses=[AIMessage("C11 is 1076.4 GPa."), AIMessage("I have no results yet; get_elastic has not run.")],
        received=[],
    )
    answer = ask(build_graph(_project(tmp_path), model), "What is C11?")
    assert answer == "I have no results yet; get_elastic has not run."
    check = model.received[1][-1]
    assert isinstance(check, HumanMessage) and "1076.4" in check.content


def test_persistent_invention_is_flagged_not_hidden(tmp_path):
    model = ScriptedModel(responses=[AIMessage("C11 is 1076.4 GPa."), AIMessage("C11 is 1076.4 GPa.")], received=[])
    answer = ask(build_graph(_project(tmp_path), model), "What is C11?")
    assert answer.startswith("C11 is 1076.4 GPa.")
    assert "[Unverified: 1076.4" in answer


def test_conversation_memory_across_turns(tmp_path):
    model = ScriptedModel(responses=[AIMessage("Hello."), AIMessage("You said hi.")], received=[])
    graph = build_graph(_project(tmp_path), model)
    ask(graph, "hi", thread_id="t")
    ask(graph, "what did I say?", thread_id="t")
    assert [m.content for m in model.received[1] if isinstance(m, HumanMessage)] == ["hi", "what did I say?"]


def test_cli_without_api_key_explains_what_to_do(tmp_path, monkeypatch, capsys):
    from elastemp.agent.cli import main

    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    monkeypatch.chdir(tmp_path)
    assert main(["--workdir", str(tmp_path), "--ask", "hi"]) == 1
    assert ".env.example" in capsys.readouterr().err


def test_cli_strips_markdown():
    from elastemp.agent.cli import plain

    text = "## Status\n**Done:** `make_elastic` ran.\n- **Next:** run VASP in `deformation-*/strain-*`."
    assert plain(text) == "Status\nDone: make_elastic ran.\n- Next: run VASP in deformation-*/strain-*."
