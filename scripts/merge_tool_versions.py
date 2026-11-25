#!/usr/bin/env python3
"""
merge_tool_versions.py - report tool table: the pipeline's tool_versions.tsv with rows replaced
by the per-sample metadata file, and the run details read from that file.

Metadata file (assets/sample_metadata/<sample>.tsv): key, value, source. Keys used by the report:
data_source, run_id, platform, library, targeting, basecalling, alignment, variants.
Keys ending in '_row' replace that tool's row in the table (value columns: version, model, path;
'basecaller' and 'aligner' rows are matched to dorado_row and minimap2_row). Same logic as the
inline step of run_panel_report.sh.

  python3 merge_tool_versions.py --tool_versions in.tsv --metadata barcode01.tsv --out out.tsv
  python3 merge_tool_versions.py --metadata barcode01.tsv --get basecalling      # one value
  python3 merge_tool_versions.py --metadata barcode01.tsv --check                # stop if incomplete
"""
import argparse
import os
import sys

REQUIRED = ("data_source", "run_id", "platform", "library", "targeting", "basecalling", "alignment", "variants")
ALIAS = {"basecaller": "dorado", "aligner": "minimap2"}


def read_meta(path):
    m = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] != "key":
            m[f[0]] = f[1:]
    return m


def problems(m):
    out = ["no value for '%s'" % k for k in REQUIRED if not (m.get(k) or [""])[0]]
    out += ["'%s' still marked %s" % (k, v[0]) for k, v in m.items()
            if v and any(t in " ".join(v).upper() for t in ("TO_CONFIRM", "UNKNOWN_SOURCE"))]
    return out


def merge(tool_rows, m):
    head, body = (tool_rows[0], tool_rows[1:]) if tool_rows else (["tool", "version", "model", "path"], [])
    over = {k[:-4]: v for k, v in m.items() if k.endswith("_row")}
    body = [r for r in body if r and ALIAS.get(r[0].lower(), r[0].lower()) not in over]
    new = [[t] + (v + ["", "", ""])[:3] for t, v in over.items()]
    return [head] + new + body


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--tool_versions", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--get", default="")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    m = read_meta(a.metadata)
    if a.check:
        p = problems(m)
        if p:
            sys.exit("[FATAL] %s: %s" % (a.metadata, "; ".join(p)))
        print("  run details complete in %s" % a.metadata)
    if a.get:
        print((m.get(a.get) or [""])[0])
    if a.out:
        rows = ([l.rstrip("\n").split("\t") for l in open(a.tool_versions)]
                if a.tool_versions and os.path.isfile(a.tool_versions) else [])
        with open(a.out, "w") as fh:
            for r in merge(rows, m):
                fh.write("\t".join(r) + "\n")
        print("  tool table: %s (rows replaced from %s: %s)" % (
              a.out, a.metadata, ", ".join(k[:-4] for k in m if k.endswith("_row")) or "none"))


if __name__ == "__main__":
    main()
