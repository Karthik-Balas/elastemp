"""Agent settings, read from the environment or a local .env file. No secrets live in the code.

Set ``ANTHROPIC_API_KEY`` (required) and optionally ``ELASTEMP_AGENT_MODEL``, either as
environment variables or in a ``.env`` file (see ``.env.example`` in the repository).
"""

import os
from pathlib import Path

API_KEY_VAR = "ANTHROPIC_API_KEY"
MODEL_VAR = "ELASTEMP_AGENT_MODEL"
DEFAULT_MODEL = "claude-sonnet-5-5"


def load_dotenv(*folders):
    """Load KEY=VALUE lines from ``.env`` in each folder into os.environ.

    Values already set in the environment are never overwritten, and empty values are ignored,
    so the empty placeholders from .env.example have no effect.
    """
    for folder in folders:
        path = Path(folder) / ".env"
        if not path.is_file():
            continue
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            key, value = key.strip(), value.strip().strip('"').strip("'")
            if value and key not in os.environ:
                os.environ[key] = value


def get_api_key():
    """The Claude API key, or a clear error explaining where to put it."""
    key = os.environ.get(API_KEY_VAR, "").strip()
    if not key:
        raise RuntimeError(
            f"No Claude API key found. Set the {API_KEY_VAR} environment variable, or copy .env.example "
            f"to .env and fill in {API_KEY_VAR}=..."
        )
    return key


def get_model(override=None):
    """Model name: command-line override, then ELASTEMP_AGENT_MODEL, then the default."""
    return override or os.environ.get(MODEL_VAR, "").strip() or DEFAULT_MODEL
