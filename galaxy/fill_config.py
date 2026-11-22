#!/usr/bin/env python3
"""
Fill a copy of the lab_strain_watch nextflow.config with user values.

The parameters of the "Parameters that must be set by the user" section of
nextflow.config are written as `name = "value" // comment` lines, inside or
outside env {} blocks. They cannot be overridden by --name on the nextflow
command line, so each web or Galaxy job gets its own filled copy of the file,
run with `nextflow -C job.config run main.nf`.

Sources of values, applied in this order (a later source wins):
    --galaxy-inputs inputs.json   the Galaxy tool form (galaxy/lab_strain_watch.xml), as
                                  written by Galaxy (<configfiles><inputs/>): sections and
                                  conditionals are nested objects, with booleans, numbers and
                                  strings. Each field named as a nextflow.config parameter
                                  sets this parameter. The other fields of the form (input
                                  type selectors, output options...) must be listed in
                                  --galaxy-ignore, so that a misspelled name is an error.
                                  The data inputs (files) are not in this JSON.
                                  --galaxy-flags FIELD=NAME1,NAME2... : FIELD is a multiple
                                  select whose options are TRUE/FALSE parameters: the
                                  selected ones are set to TRUE, the other NAMEs to FALSE
                                  (not used by the lab_strain_watch form).
    --values params.json          JSON object {name: value} typed by a user (e.g. the web page
                                  or the optional parameter file of the Galaxy tool).
    --set name=value              values set by the server or the Galaxy wrapper (trusted,
                                  e.g. "$PWD/results").

Values of the first two sources are untrusted and strictly checked: no quote,
backslash, dollar, backtick or newline, so that nothing can be interpolated or
executed by Groovy when nextflow reads the config. Booleans become "TRUE" or
"FALSE" and numbers are written without exponent ("1" for 1.0), as expected by
conf/CheckVariables.nf.

--cpus and --memory-mb add the resources of the Galaxy job at the end of the
config (executor.cpus/memory of the nextflow local executor and
process.resourceLimits, which lowers the cpus/memory of the process labels).

Usage:
    fill_config.py --template nextflow.config --spec > spec.json
    fill_config.py --template nextflow.config --galaxy-inputs inputs.json \
        --galaxy-ignore input_type,sample_name --values params.json \
        --set system_exec=local --set out_path_ini=/tmp/res --cpus 8 -o job.config
"""

import argparse
import decimal
import json
import math
import re
import sys

PARAM_RE = re.compile(
    r"^(?P<indent>\s*)(?P<name>[A-Za-z_][A-Za-z0-9_]*)(?P<eq>\s*=\s*)"
    r"(?P<value>\"[^\"]*\"|'[^']*'|[0-9]+)(?P<rest>.*)$"
)
SECTION_RE = re.compile(r"^##\s+(?P<title>\S.*?)\s*##$")
END_USER_SECTION = "End Parameters that must be set by the user"
FORBIDDEN_CHARS = set("\"'\\$\n\r`")


def parse(lines):
    """Return the list of user parameters of nextflow.config, in file order.

    Each parameter is a dict with: name, default, quote, section, help, line.
    """
    params = []
    section = ""
    in_comment_block = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if END_USER_SECTION in stripped:
            break
        if stripped.startswith("/*"):
            in_comment_block = True
        if in_comment_block:
            m = SECTION_RE.match(stripped)
            if m and "Parameters that must be set" not in m.group("title"):
                section = m.group("title")
            if stripped.endswith("*/"):
                in_comment_block = False
            continue
        m = PARAM_RE.match(line.rstrip("\r\n"))
        if m:
            value = m.group("value")
            quote = value[0] if value[0] in "\"'" else ""
            rest = m.group("rest").strip()
            help_text = rest[2:].strip() if rest.startswith("//") else ""
            params.append({
                "name": m.group("name"),
                "default": value[1:-1] if quote else value,
                "quote": quote,
                "section": section,
                "help": help_text,
                "line": i,
            })
        elif stripped.startswith("//") and params and params[-1]["line"] < i and not stripped.startswith("///"):
            # continuation lines of the previous parameter comment, only inside env {} (indented lines right after the parameter)
            prev = lines[params[-1]["line"] + 1:i]
            if all(p.strip().startswith("//") for p in prev) and line.startswith((" ", "\t")):
                params[-1]["help"] += "\n" + stripped[2:].strip()
    return params


def check_user_value(name, value):
    if not isinstance(value, (str, int, float)) or isinstance(value, bool):
        raise ValueError(f"{name}: value must be a string or a number")
    value = str(value)
    bad = FORBIDDEN_CHARS.intersection(value)
    if bad:
        raise ValueError(f"{name}: forbidden character(s) {''.join(sorted(bad))!r}")
    if len(value) > 1000:
        raise ValueError(f"{name}: value too long")
    return value


def format_number(name, value):
    """Number of the Galaxy form -> string accepted by conf/CheckVariables.nf (no exponent, 1 and not 1.0)."""
    if isinstance(value, int):
        return str(value)
    if not math.isfinite(value):
        raise ValueError(f"{name}: invalid number {value!r}")
    if value.is_integer():
        return str(int(value))
    text = format(decimal.Decimal(repr(value)), "f")
    return text.rstrip("0").rstrip(".") if "." in text else text


def galaxy_value(name, value):
    """Value of the Galaxy inputs JSON -> string of nextflow.config."""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, (int, float)):
        return format_number(name, value)
    if isinstance(value, str):
        return value
    raise ValueError(f"{name}: unexpected value {value!r} in the Galaxy inputs (list it in --galaxy-ignore if it is not a nextflow.config parameter)")


def flatten_galaxy_inputs(tree, ignore, where=""):
    """Return {name: value} of the fields of the Galaxy inputs JSON. Sections and conditionals are nested objects."""
    out = {}
    for key, value in tree.items():
        path = f"{where}|{key}" if where else key
        if key in ignore:
            continue
        if isinstance(value, dict):
            items = flatten_galaxy_inputs(value, ignore, path).items()
        else:
            items = [(key, (path, value))]
        for name, item in items:
            if name in out:
                raise ValueError(f"{name}: field present twice in the Galaxy form ({out[name][0]} and {item[0]})")
            out[name] = item
    return out if where else {name: value for name, (path, value) in out.items()}


def fill(lines, params, values):
    """Return the lines with the values replaced. values: {name: str}"""
    by_name = {p["name"]: p for p in params}
    unknown = set(values) - set(by_name)
    if unknown:
        raise ValueError("unknown parameter(s): " + ", ".join(sorted(unknown)))
    out = list(lines)
    for name, value in values.items():
        p = by_name[name]
        if p["quote"]:
            if p["quote"] in value or "\n" in value:
                raise ValueError(f"{name}: value cannot contain {p['quote']} or a newline")
            new_value = p["quote"] + value + p["quote"]
        elif re.fullmatch(r"[0-9]+", value):
            new_value = value
        else:
            raise ValueError(f"{name}: integer expected")
        line = out[p["line"]]
        eol = line[len(line.rstrip("\r\n")):]
        m = PARAM_RE.match(line.rstrip("\r\n"))
        out[p["line"]] = m.group("indent") + name + m.group("eq") + new_value + m.group("rest") + eol
    return out


def resource_lines(cpus, memory_mb):
    """Lines appended to the config: resources of the Galaxy job for the nextflow local executor."""
    if not cpus and not memory_mb:
        return ""
    lines = ["", "", "// Added by galaxy/fill_config.py: resources of the Galaxy job, used by the local executor of nextflow"]
    limits = []
    if cpus:
        lines.append(f"executor.cpus = {cpus}")
        limits.append(f"cpus: {cpus}")
    if memory_mb:
        lines.append(f"executor.memory = '{memory_mb} MB'")
        limits.append(f"memory: {memory_mb}.MB")
    if limits:
        lines.append(f"process.resourceLimits = [{', '.join(limits)}]  // lowers the cpus/memory of the process labels if above")
    return "\n".join(lines) + "\n"


def positive_int_or_empty(text):
    """argparse type: "" (unset variable in the Galaxy command) -> None."""
    text = text.strip()
    if text == "":
        return None
    if not re.fullmatch(r"[0-9]+", text) or int(text) < 1:
        raise argparse.ArgumentTypeError(f"positive integer expected, not {text!r}")
    return int(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--template", required=True, help="path of the nextflow.config of lab_strain_watch")
    parser.add_argument("--spec", action="store_true", help="print the user parameters as JSON and exit")
    parser.add_argument("--galaxy-inputs", help="JSON file of the Galaxy tool form (<configfiles><inputs/> of the tool)")
    parser.add_argument("--galaxy-ignore", action="append", default=[], metavar="NAME[,NAME...]",
                        help="fields of the Galaxy form that are not nextflow.config parameters, can be repeated")
    parser.add_argument("--galaxy-flags", action="append", default=[], metavar="FIELD=NAME1,NAME2...",
                        help="multiple select FIELD of the Galaxy form whose options are TRUE/FALSE parameters, can be repeated")
    parser.add_argument("--values", help="JSON file of untrusted user values {name: value}")
    parser.add_argument("--set", action="append", default=[], metavar="NAME=VALUE", help="trusted value, can be repeated")
    parser.add_argument("--locked", help="file of parameter names that users cannot set (one per line, # for comments)")
    parser.add_argument("--cpus", type=positive_int_or_empty, help="cpus of the job (empty: not set)")
    parser.add_argument("--memory-mb", type=positive_int_or_empty, help="memory of the job in MB (empty: not set)")
    parser.add_argument("-o", "--output", help="output file (default: stdout)")
    args = parser.parse_args()

    with open(args.template, encoding="utf-8", newline="") as f:
        lines = f.readlines()
    params = parse(lines)

    if args.spec:
        json.dump(params, sys.stdout, indent=1)
        sys.stdout.write("\n")
        return

    locked = set()
    if args.locked:
        with open(args.locked, encoding="utf-8") as f:
            locked = {l.split("#")[0].strip() for l in f} - {""}

    values = {}
    try:
        if args.galaxy_inputs:
            with open(args.galaxy_inputs, encoding="utf-8") as f:
                tree = json.load(f)
            if not isinstance(tree, dict):
                raise ValueError("the --galaxy-inputs file must contain a JSON object")
            ignore = {n.strip() for item in args.galaxy_ignore for n in item.split(",")} - {""}
            form = flatten_galaxy_inputs(tree, ignore)
            for item in args.galaxy_flags:
                field, sep, names = item.partition("=")
                names = [n.strip() for n in names.split(",") if n.strip()]
                if not sep or not names:
                    raise ValueError(f"--galaxy-flags {item}: FIELD=NAME1,NAME2... expected")
                if field not in form:  # field of a conditional case not selected
                    continue
                selected = form.pop(field) or []
                if not isinstance(selected, list) or not set(selected) <= set(names):
                    raise ValueError(f"{field}: unexpected value {selected!r} in the Galaxy inputs")
                for name in names:
                    if name in form:
                        raise ValueError(f"{name}: set by the {field} field and by another field of the Galaxy form")
                    form[name] = "TRUE" if name in selected else "FALSE"
            forbidden = locked.intersection(form)
            if forbidden:
                raise ValueError("locked parameter(s) in the Galaxy form: " + ", ".join(sorted(forbidden)))
            for name, value in form.items():
                if value is None:  # optional field left empty: nextflow.config default
                    continue
                values[name] = check_user_value(name, galaxy_value(name, value))
        if args.values:
            with open(args.values, encoding="utf-8") as f:
                user_values = json.load(f)
            if not isinstance(user_values, dict):
                raise ValueError("the --values file must contain a JSON object")
            forbidden = locked.intersection(user_values)
            if forbidden:
                raise ValueError("locked parameter(s): " + ", ".join(sorted(forbidden)))
            values.update({k: check_user_value(k, v) for k, v in user_values.items()})
        for item in args.set:
            name, sep, value = item.partition("=")
            if not sep:
                raise ValueError(f"--set {item}: NAME=VALUE expected")
            values[name] = value
        out = fill(lines, params, values)
    except (ValueError, json.JSONDecodeError) as e:
        sys.exit(f"fill_config.py: ERROR: {e}")

    text = "".join(out) + resource_lines(args.cpus, args.memory_mb)
    if args.output:
        with open(args.output, "w", encoding="utf-8", newline="") as f:
            f.write(text)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
