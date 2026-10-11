import re
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import pytest

import havelock as hv
from havelock.nonlinear import c4_level_set, nonlinear

REPO = Path(__file__).resolve().parents[1]

# paper Table c4 (level-set column), n = 4..25
C4_TABLE = {4: "1/4", 5: "9/40", 6: "23/30", 7: "9/14", 8: "95/56", 9: "19/12", 10: "49/15", 11: "75/22",
            12: "755/132", 13: "675/104", 14: "1693/182", 15: "56/5", 16: "1141/80", 17: "609/34",
            18: "3190/153", 19: "27", 20: "11139/380", 21: "2175/56", 22: "6145/154", 23: "2475/46",
            24: "29183/552", 25: "7227/100"}


@pytest.mark.parametrize("n", sorted(C4_TABLE))
def test_c4_table(n):
    r = nonlinear("mt", n)
    assert r.c4 == Fraction(C4_TABLE[n]) and r.verdict == "stable" and r.status == hv.Status.PROVED
    assert r.mode == n // 2 and r.method == "closed form" and "Havelock.Conj13" in r.basis


def test_views():
    assert nonlinear("plane", 7).c4 == Fraction(9, 14)
    assert nonlinear("sinh", 8).c4 == Fraction(95, 56)
    assert nonlinear("sphere", 4).c4 == Fraction(1, 4)


@pytest.mark.parametrize("model,N", [("mt", 3), ("sphere", 9), ("plane", 8), ("sinh", 5), ("mt", 2)])
def test_out_of_validity(model, N):
    with pytest.raises(hv.OutOfValidityRange):
        nonlinear(model, N)


def test_reproduces_recorded_table_check():
    # research reproduction (spec 9.7): conj13/paper/table_check.out records the closed-form c4 for each n
    text = (REPO / "conj13/paper/table_check.out").read_text()
    rows = re.findall(r"^n=(\d+) b_n=\S+ c4\(float\)=\S+ formula=(\d+(?:/\d+)?)=", text, re.M)
    assert len(rows) >= 15
    for n, frac in rows:
        assert nonlinear("mt", int(n)).c4 == Fraction(frac)


def test_level_set_engine_reproduces_closed_forms():
    # verified while planning: |engine - 95/56| = 1.2e-24 (n = 8), |engine - 1/4| = 1.2e-26 (n = 4) at dps 40
    with mp.workdps(50):
        b8 = -(8 - 3 * mp.sqrt(7))
        b4 = 2 - mp.sqrt(3)
        assert abs(c4_level_set("mt", 8, b8) - mp.mpf(95) / 56) < mp.mpf("1e-18")
        assert abs(c4_level_set("mt", 4, b4) - mp.mpf(1) / 4) < mp.mpf("1e-18")


def test_level_set_engine_odd_n():
    # n = 7 (plane): verified while planning, |engine - 9/14| = 3.3e-17 at dps 50
    with mp.workdps(50):
        assert abs(c4_level_set("mt", 7, 0, dps=50) - mp.mpf(9) / 14) < mp.mpf("1e-12")
