"""``elastemp-agent``: chat with the elastemp agent about one project folder.

    elastemp-agent --workdir runs/TiN                      # interactive
    elastemp-agent --workdir runs/TiN --ask "What's next?"  # one question

Needs ANTHROPIC_API_KEY in the environment or in a .env file (see .env.example).
"""

import argparse
import logging
import re
import sys
import uuid

from elastemp.base.utils import resolve_workspace


def plain(text):
    """Strip Markdown emphasis, inline code and headings in case the model uses them anyway."""
    text = re.sub(r"\*\*(.+?)\*\*", r"\1", text)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    return re.sub(r"^#{1,6}\s+", "", text, flags=re.MULTILINE)


def _configure_file_logging(workspace):
    """elastemp's own log goes to the project's elastemp.log only, keeping the chat readable."""
    logger = logging.getLogger("elastemp")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    handler = logging.FileHandler(workspace / "elastemp.log")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    logger.addHandler(handler)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="elastemp-agent", description="Chat with the elastemp agent.")
    parser.add_argument("--workdir", default=".", help="project folder (default: current directory)")
    parser.add_argument("--model", default=None, help="Claude model (default: ELASTEMP_AGENT_MODEL or built-in)")
    parser.add_argument("--ask", default=None, help="ask one question and exit")
    parser.add_argument("--quiet", action="store_true", help="do not show tool calls")
    args = parser.parse_args(argv)

    workspace = resolve_workspace(args.workdir)
    if not workspace.is_dir():
        parser.error(f"--workdir {workspace} is not a directory")

    try:
        from elastemp.agent.config import load_dotenv
        from elastemp.agent.graph import ask, build_graph, make_claude
    except ImportError as exc:
        print(f'The agent needs extra packages: pip install "elastemp[agent]"  ({exc})', file=sys.stderr)
        return 1

    load_dotenv(".", workspace)
    try:
        llm = make_claude(args.model)
    except RuntimeError as exc:
        print(exc, file=sys.stderr)
        return 1

    _configure_file_logging(workspace)
    graph = build_graph(workspace, llm)
    thread = str(uuid.uuid4())

    def show(kind, detail):
        if not args.quiet:
            print(f"  [{kind}] {detail}", file=sys.stderr)

    if args.ask:
        print(plain(ask(graph, args.ask, thread, show)))
        return 0

    print(f"elastemp agent for {workspace}. Type 'exit' to quit.")
    while True:
        try:
            question = input("\nyou> ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return 0
        if question.lower() in {"exit", "quit"}:
            return 0
        if question:
            print(f"\nagent> {plain(ask(graph, question, thread, show))}")


if __name__ == "__main__":
    sys.exit(main())
