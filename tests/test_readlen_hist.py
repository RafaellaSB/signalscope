#!/usr/bin/env python3
"""Tests for scripts/readlen_hist.py (python3 tests/test_readlen_hist.py)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import readlen_hist as rh  # noqa: E402


def test_bins_are_log_spaced():
    assert rh.bin_of(1000) == 60 and rh.bin_of(10000) == 80 and rh.bin_of(999) == 59


def test_histogram_from_samtools_stats_rl():
    lines = ["SN\traw total sequences:\t5", "RL\t1000\t3", "RL\t10000\t2", "RL\t0\t9"]
    assert dict(rh.histogram(lines)) == {60: 3, 80: 2}


if __name__ == "__main__":
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for k, f in tests:
        f()
        print("ok", k)
    print("%d tests passed" % len(tests))
