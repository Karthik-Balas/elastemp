"""The elastemp agent: a LangGraph ReAct loop with Claude, plus a number-grounding check.

Graph:

    START -> agent --(tool calls)--> tools -> agent
                   --(final answer)--> verify --(ok)--> END
                                              --(invented numbers, first time)--> agent
                                              --(still invented)--> END, with a warning attached

The agent never runs VASP and never submits jobs; the user does. Steps that would delete
finished VASP runs are refused by the run_step tool unless the user explicitly confirms.
"""

import json

from langchain_core.messages import AIMessage, HumanMessage, SystemMessage, ToolMessage
from langchain_core.tools import tool
from langgraph.checkpoint.memory import InMemorySaver
from langgraph.graph import END, START, MessagesState, StateGraph
from langgraph.prebuilt import ToolNode

from elastemp.agent import tools as T
from elastemp.agent.grounding import ungrounded_numbers
from elastemp.base.utils import resolve_workspace

MAX_VERIFY_RETRIES = 1
CHECK_PREFIX = "[automatic number check]"

SYSTEM_PROMPT = """You are the assistant for elastemp, a workflow that computes temperature-dependent elastic \
constants of a crystal from VASP and phonopy calculations, using the quasi-harmonic approximation.

The project folder is: {workspace}

Workflow, in order (the user runs VASP between steps; you never run or submit VASP jobs):
1. make_elastic: create strained structures. User runs VASP in deformation-*/strain-*.
2. get_elastic: fit energy vs strain; 0 K elastic constants, moduli, stability.
3. make_bulk_dynamic: phonon supercells for deformation-bulk. User runs VASP DFPT there.
4. get_bulk_dynamic: thermal expansion from the bulk phonons.
5. make_temp_static: structures scaled to the tmax volume. User runs VASP in Hightemp_static.
6. get_temp_static: energies at the tmax volume.
7. make_dynamic: phonon supercells for every deformation. User runs VASP DFPT there.
8. get_dynamic: force constants.
9. get_thermal: elastic constants and moduli from 0 K to tmax.

Your tools:
- project_status: which steps are done and what comes next. Call it before answering any question about \
progress or what to do next.
- run_step: run one step. Only run a step when the user asks you to, or agrees to your suggestion. If it \
returns status "refused", explain what would be deleted and ask the user; never set confirm=true unless the \
user has explicitly confirmed in this conversation.
- inspect_run: evidence for one VASP run (convergence, OSZICAR, OUTCAR warnings). Use it to explain failures.
- read_results: all computed numbers and consistency checks.

Rules:
- Every number you state must come from a tool result in this conversation. Do not calculate new numbers, \
and do not quote values from memory or the literature (for example experimental elastic constants). If a \
number you need is not in a tool result, say so.
- Always mention every item in "flags" from read_results, and any failed or dropped runs.
- Be concise and concrete. Use plain language a materials scientist would use; give units.
- Your answer is shown in a plain terminal: write plain text. No Markdown (no **bold**, no backticks, \
no # headings); use simple "-" bullets if you need a list.
"""


def _text(message):
    content = message.content
    if isinstance(content, str):
        return content
    return "".join(block.get("text", "") for block in content if isinstance(block, dict))


class AgentState(MessagesState):
    verify_attempts: int


def make_tools(workspace):
    """The four tools, bound to one project folder."""
    ws = str(resolve_workspace(workspace))

    @tool
    def project_status() -> str:
        """Report which workflow steps are done, which VASP runs are finished, failed or missing,
        and the next action."""
        return T.tool_project_status(ws)

    @tool
    def run_step(step: str, confirm: bool = False) -> str:
        """Run one elastemp workflow step in the project folder.

        Args:
            step: one of make_elastic, get_elastic, make_bulk_dynamic, get_bulk_dynamic, make_temp_static,
                get_temp_static, make_dynamic, get_dynamic, get_thermal.
            confirm: set to true only after the user has explicitly agreed to a step that would delete
                existing VASP results.
        """
        return T.tool_run_step(ws, step, confirm)

    @tool
    def inspect_run(folder: str) -> str:
        """Show convergence status and OSZICAR/OUTCAR evidence for one VASP run folder.

        Args:
            folder: path relative to the project folder, e.g. deformation-2/strain-4 or
                deformation-bulk/strain-3/Phonon_calculation.
        """
        return T.tool_inspect_run(ws, folder)

    @tool
    def read_results() -> str:
        """Return all computed results (elastic constants, moduli, stability, temperature dependence),
        consistency checks and flags."""
        return T.tool_read_results(ws)

    return [project_status, run_step, inspect_run, read_results]


def build_graph(workspace, llm, *, checkpointer=None):
    """Compile the agent graph.

    :param llm: a LangChain chat model that supports tool calling (ChatAnthropic in normal use)
    :param checkpointer: LangGraph checkpointer for multi-turn memory (default: in-memory)
    """
    ws = resolve_workspace(workspace)
    tools = make_tools(ws)
    llm_with_tools = llm.bind_tools(tools)
    system = SystemMessage(SYSTEM_PROMPT.format(workspace=ws))

    def agent(state: AgentState):
        return {"messages": [llm_with_tools.invoke([system, *state["messages"]])]}

    def route_after_agent(state: AgentState):
        return "tools" if state["messages"][-1].tool_calls else "verify"

    def verify(state: AgentState):
        answer = state["messages"][-1]
        references = [
            _text(m) for m in state["messages"] if isinstance(m, ToolMessage | HumanMessage) and not _is_check(m)
        ]
        bad = ungrounded_numbers(_text(answer), references)
        if not bad:
            return {}
        attempts = state.get("verify_attempts", 0)
        if attempts < MAX_VERIFY_RETRIES:
            note = (
                f"{CHECK_PREFIX} These numbers in your answer do not appear in any tool result: {', '.join(bad)}. "
                "Rewrite the answer using only values from tool results (call a tool if you need one), or remove them."
            )
            return {"messages": [HumanMessage(note)], "verify_attempts": attempts + 1}
        warning = f"\n\n[Unverified: {', '.join(bad)} could not be traced to a tool result.]"
        return {"messages": [AIMessage(content=_text(answer) + warning, id=answer.id)]}

    def route_after_verify(state: AgentState):
        last = state["messages"][-1]
        return "agent" if isinstance(last, HumanMessage) and _is_check(last) else END

    graph = StateGraph(AgentState)
    graph.add_node("agent", agent)
    graph.add_node("tools", ToolNode(tools))
    graph.add_node("verify", verify)
    graph.add_edge(START, "agent")
    graph.add_conditional_edges("agent", route_after_agent, ["tools", "verify"])
    graph.add_edge("tools", "agent")
    graph.add_conditional_edges("verify", route_after_verify, ["agent", END])
    return graph.compile(checkpointer=checkpointer or InMemorySaver())


def _is_check(message):
    return isinstance(message, HumanMessage) and _text(message).startswith(CHECK_PREFIX)


def make_claude(model=None, max_tokens=4096):
    """ChatAnthropic configured from the environment / .env (see elastemp.agent.config)."""
    from langchain_anthropic import ChatAnthropic

    from elastemp.agent.config import get_api_key, get_model

    return ChatAnthropic(model=get_model(model), api_key=get_api_key(), max_tokens=max_tokens)


def ask(graph, question, thread_id="default", on_event=None):
    """Send one user message; return the final answer text. ``on_event`` receives (kind, detail) updates."""
    config = {"configurable": {"thread_id": thread_id}, "recursion_limit": 40}
    inputs = {"messages": [HumanMessage(question)], "verify_attempts": 0}
    for update in graph.stream(inputs, config, stream_mode="updates"):
        if on_event is None:
            continue
        for node, payload in update.items():
            for message in (payload or {}).get("messages", []):
                if node == "agent" and getattr(message, "tool_calls", None):
                    for call in message.tool_calls:
                        args = json.dumps(call["args"]) if call["args"] else ""
                        on_event("tool", f"{call['name']}({args})")
                elif node == "verify" and _is_check(message):
                    on_event("check", _text(message))
    return _text(graph.get_state(config).values["messages"][-1])
