import os
import re
import json
from typing import Any

# regex to capture the floats in your naming scheme
# e.g. FittingTest_s[0.0000_2.0000]_T[-100.0000_-200.0000]_a[1.0000]_noise[0.0000]_Rot[0.0000].png
FILENAME_PATTERN = re.compile(
    r"""FittingTest
        _s\[
          (?P<s0>[-\d\.]+)_
          (?P<s1>[-\d\.]+)
        \]
        _T\[
          (?P<T0>[-\d\.]+)_
          (?P<T1>[-\d\.]+)
        \]
        _a\[(?P<a>[-\d\.]+)\]
        _noise\[(?P<noise>[-\d\.]+)\]
        _Rot\[(?P<rot>[-\d\.]+)\]
        \.png$
    """,
    re.VERBOSE,
)


def parse_filename(
    fn: str,
) -> tuple[tuple[float, float], tuple[float, float], float, float, float]:
    """Extract (s0,s1), (T0,T1), a, noise, rot from a filename."""
    m = FILENAME_PATTERN.match(fn)
    if not m:
        raise ValueError(f"Filename {fn!r} doesn't match pattern")
    s0, s1 = float(m.group("s0")), float(m.group("s1"))
    T0, T1 = float(m.group("T0")), float(m.group("T1"))
    a = float(m.group("a"))
    noise = float(m.group("noise"))
    rot = float(m.group("rot"))
    return (s0, s1), (T0, T1), a, noise, rot


def generate_index(img_dir: str, out_json: str) -> None:
    """
    Scan img_dir for PNGs, group by (s,T,a,noise),
    sort by rot, and write out a JSON list of groups.
    """
    print(f"Generating index.json from {img_dir} to {out_json}", end="")

    groups: dict[tuple, list[dict[str, Any]]] = {}
    for fn in os.listdir(img_dir):
        if not fn.endswith(".png"):
            continue
        key = parse_filename(fn)[:4]  # (s,T,a,noise)
        rot = parse_filename(fn)[4]
        groups.setdefault(key, []).append({"filename": fn, "rot": rot})

    # sort each list by rot
    out_list: list[dict[str, Any]] = []
    for (s, T, a, noise), imgs in groups.items():
        imgs.sort(key=lambda x: x["rot"])
        out_list.append({"s": s, "T": T, "a": a, "noise": noise, "images": imgs})

    # write JSON
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(out_list, f, indent=2)

    print("Complete")


if __name__ == "__main__":
    # adjust paths as needed
    generate_index(img_dir="img", out_json="index.json")
