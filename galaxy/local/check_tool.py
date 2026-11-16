#!/usr/bin/env python3
"""
Check that the form of the Galaxy tool (galaxy/biosentinelle.xml) still recapitulates nextflow.config.

Errors (exit code 1):
    - a user parameter of nextflow.config (not locked) is not in the tool (form field, option value or command)
    - a field of the form is neither a nextflow.config parameter, nor a data input, nor in the galaxy-ignore list of the command
    - the default of a select is not one of its options
Warnings:
    - the default of a form field differs from the nextflow.config default (numbers compared as numbers: 0.9 is 0.90)

Run by galaxy/local/test.sh, after any change of nextflow.config or of the tool.

Usage:
    check_tool.py --template nextflow.config --locked galaxy/locked_params.txt --tool galaxy/biosentinelle.xml
"""

import argparse
import decimal
import os
import re
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fill_config import format_number, parse  # noqa: E402


def form_default(param):
    """Default of a form field, as written in nextflow.config (None if not comparable)."""
    kind = param.get("type")
    if kind == "boolean":
        return "TRUE" if param.get("checked", "false").lower() in ("true", "yes") else "FALSE"
    if kind == "select":
        options = param.findall("option")
        selected = [o.get("value") for o in options if o.get("selected", "false").lower() in ("true", "yes")]
        return selected[0] if selected else (options[0].get("value") if options else None)
    if kind in ("integer", "float"):
        value = param.get("value")
        return None if value in (None, "") else format_number(param.get("name"), float(value) if kind == "float" else int(value))
    if kind in ("text", "hidden"):
        return param.get("value", "")
    return None


def same_value(a, b):
    """Same default: equal strings, or equal numbers (0.9 and 0.90)."""
    if a == b:
        return True
    try:
        return decimal.Decimal(a) == decimal.Decimal(b)
    except (decimal.InvalidOperation, TypeError):
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--template", required=True, help="nextflow.config of biosentinelle")
    parser.add_argument("--locked", required=True, help="parameters not in the form (galaxy/locked_params.txt)")
    parser.add_argument("--tool", required=True, help="galaxy/biosentinelle.xml")
    args = parser.parse_args()

    with open(args.template, encoding="utf-8") as f:
        config = {p["name"]: p for p in parse(f.readlines())}
    with open(args.locked, encoding="utf-8") as f:
        locked = {l.split("#")[0].strip() for l in f} - {""}
    with open(args.tool, encoding="utf-8") as f:
        tool_text = f.read()
    tool = ET.fromstring(tool_text)
    command = tool.findtext("command") or ""
    m = re.search(r"--galaxy-ignore\s+(\S+)", command)
    ignore = set(m.group(1).split(",")) if m else set()
    # multiple selects of TRUE/FALSE parameters: --galaxy-flags FIELD=NAME1,NAME2...
    ignore |= {field for field in re.findall(r"--galaxy-flags\s+([A-Za-z0-9_]+)=", command)}

    errors, warnings = [], []
    params = tool.find("inputs").iter("param")
    form = {}
    for param in params:
        name = param.get("name")
        form.setdefault(name, []).append(param)
        if name in config:
            if name in locked:
                errors.append(f"{name}: locked parameter in the form")
            continue
        if param.get("type") in ("data", "data_collection") or name in ignore:
            continue
        errors.append(f"{name}: form field that is not a nextflow.config parameter (add it to the galaxy-ignore option of the command if it is a Galaxy-only field)")

    for name, items in form.items():
        for param in items:
            if param.get("type") == "select":
                values = [o.get("value") for o in param.findall("option")]
                default = form_default(param)
                if default is not None and default not in values:
                    errors.append(f"{name}: default {default!r} is not an option")
            if name in config:
                default = form_default(param)
                if default is not None and not same_value(default, config[name]["default"]):
                    warnings.append(f"{name}: default of the form {default!r} differs from nextflow.config {config[name]['default']!r}")

    for name, p in config.items():
        if name in locked or name in form:
            continue
        # parameters set by the command (e.g., options of a select)
        if not re.search(r"\b" + re.escape(name) + r"\b", tool_text):
            errors.append(f"{name}: nextflow.config parameter ({p['section']} section) missing from the tool and not locked")

    for name in sorted(locked - set(config)):
        warnings.append(f"{name}: in {os.path.basename(args.locked)} but not a nextflow.config parameter")

    for w in warnings:
        print(f"WARNING: {w}")
    for e in errors:
        print(f"ERROR: {e}")
    n_form = len([n for n in config if n in form])
    print(f"{len(config)} parameters in nextflow.config: {n_form} in the form, {len(locked & set(config))} locked, "
          f"{len(config) - n_form - len(locked & set(config))} set by the command. {len(errors)} error(s), {len(warnings)} warning(s).")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
