#!/usr/bin/env python3
"""Tests for the long-read SV frequency parser and matching in sv_dosage_annotate.py."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import sv_dosage_annotate as sd  # noqa: E402

# two CoLoRSdb-style records (sequence-resolved, AF only in the record) and one with AC/AN only
LINES = ["chr1\t11028\t.\t.\tSVTYPE=DEL;SVLEN=-35;AC=50;AN=2854;NS=1427;AF=0.0175193;nhomalt=1",
         "chr1\t20000\t.\t.\tSVTYPE=INS;SVLEN=300;AC=10;AN=2854",
         "chr2\t500\t.\t.\tSVTYPE=DEL;SVLEN=-1000;AC=2854;AN=2854"]


def test_parse_af_end_and_acan():
    gn = sd.parse_lr_records(LINES)
    d = list(gn[("1", "DEL")].values())[0]
    assert d[0] == 11028 and d[1] == 11028 + 35 and abs(d[3] - 0.0175193) < 1e-9
    i = list(gn[("1", "INS")].values())[0]
    assert i[1] == 20000 and i[5] == 300 and abs(i[3] - 10 / 2854) < 1e-9
    assert list(gn[("2", "DEL")].values())[0][3] == 1.0


def test_match_del_by_overlap_and_ins_by_distance_and_size():
    gn = sd.parse_lr_records(LINES)
    svs = [{"chrom": "1", "type": "DEL", "start": 11030, "end": 11064, "len": 34},
           {"chrom": "1", "type": "INS", "start": 20300, "end": 20300, "len": 260},
           {"chrom": "1", "type": "INS", "start": 20300, "end": 20300, "len": 100}]
    m = sd.match_gnomad(list(enumerate(svs)), gn, 0.5, 500, 0.7)
    assert 0 in m and 1 in m and 2 not in m          # size ratio 100/300 < 0.7


def test_malformed_records_skipped():
    gn = sd.parse_lr_records(["chr1\t5\t.\t.\tSVTYPE=DEL;SVLEN=x", "bad line"])
    assert not gn


if __name__ == "__main__":                      # runs without pytest too
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for k, f in tests:
        f()
        print("ok", k)
    print("%d tests passed" % len(tests))
