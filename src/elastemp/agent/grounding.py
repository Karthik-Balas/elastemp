"""Check that every number in the agent's answer comes from a tool result (or the user).

The model is told to quote numbers only from tools. This module checks that it did, so an
invented value (say "C11 = 1076 GPa") is caught by code rather than trusted.
"""

import re

# Identifiers that contain digits but are not quantities: C11, c44, strain-4, deformation-3, NELM...
_IDENTIFIER = re.compile(r"\b[A-Za-z][A-Za-z_]*[-_]?\d+\b")
_NUMBER = re.compile(r"(?<![\w.])[-−]?\d+(?:,\d{3})*(?:\.\d+)?(?:[eE][-+]?\d+)?(?![\w])")
_PLAIN_NUMBER = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")

SMALL_INTEGER_LIMIT = 30  # counts, step numbers, list markers: not worth checking


def numbers_in_answer(text):
    """Numbers stated in prose, as (text, value) pairs, skipping identifiers and small integers."""
    cleaned = _IDENTIFIER.sub(" ", text)
    found = []
    for match in _NUMBER.finditer(cleaned):
        raw = match.group(0)
        value = float(raw.replace(",", "").replace("−", "-"))
        if "." not in raw and "e" not in raw.lower() and abs(value) <= SMALL_INTEGER_LIMIT:
            continue
        found.append((raw, value))
    return found


def reference_numbers(texts):
    """Every number appearing in the given texts (tool outputs, user messages)."""
    values = set()
    for text in texts:
        for raw in _PLAIN_NUMBER.findall(text):
            try:
                values.add(float(raw))
            except ValueError:
                pass
    return values


def _decimals(raw):
    raw = raw.lower().split("e")[0]
    return len(raw.split(".")[1]) if "." in raw else 0


def is_grounded(raw, value, references):
    """True if ``value`` matches a reference once rounding to the stated precision is allowed.

    The sign is ignored, so "dropped by 4.5%" matches a stored change of -4.54.
    """
    tolerance = 0.5 * 10 ** (-_decimals(raw))
    target = abs(value)
    return any(abs(target - abs(ref)) <= max(tolerance, 1e-12) * 1.0001 for ref in references)


def ungrounded_numbers(answer, reference_texts):
    """Numbers in ``answer`` that do not appear, up to rounding, in any reference text."""
    references = reference_numbers(reference_texts)
    return [raw for raw, value in numbers_in_answer(answer) if not is_grounded(raw, value, references)]
