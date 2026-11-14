import re
import sys
from pathlib import Path

MODELS = Path(__file__).resolve().parent / "models"
# (arquivo, rotulo, alvo L, alvo W) - bin escolhido para as celulas em cells/
TARGETS = [
    ("nfet_tt.pm3.spice", "NFET", 0.150e-6, 1.0e-6),
    ("pfet_tt.pm3.spice", "PFET", 0.150e-6, 1.0e-6),
]


def find_bin(text: str, l: float, w: float):
    """Return the first lmin/lmax/wmin/wmax bin that contains (l, w)."""
    bins = re.finditer(
        r"\+ lmin = ([\d.e+\-]+) lmax = ([\d.e+\-]+) "
        r"wmin = ([\d.e+\-]+) wmax = ([\d.e+\-]+)",
        text,
    )
    for m in bins:
        lmin, lmax, wmin, wmax = map(float, m.groups())
        if lmin <= l <= lmax and wmin <= w <= wmax:
            return m, lmin, lmax, wmin, wmax
    return None


PARAMS = (
    ("vth0", r"\+ vth0 = \{?\{?([\d.e\-]+)"),
    ("u0", r"\+ u0 = ([\d.e\-]+)"),
    ("k1", r"\+ k1 = ([\d.e\-]+)"),
    ("ndep", r"\+ ndep = ([\d.e+]+17)"),
    ("toxm", r"\+ toxm = ([\d.e\-]+)"),
)


def main() -> int:
    for fname, label, l, w in TARGETS:
        text = (MODELS / fname).read_text(encoding="utf-8", errors="ignore")
        hit = find_bin(text, l, w)
        if hit is None:
            print(f"{label}: no bin for L={l:.2e} W={w:.2e}")
            continue
        m, lmin, lmax, wmin, wmax = hit
        block = text[m.start() : m.start() + 60000]
        print(f"{label}: bin L[{lmin:.2e},{lmax:.2e}] W[{wmin:.2e},{wmax:.2e}]")
        for name, pattern in PARAMS:
            found = re.search(pattern, block)
            print(f"  {name} = {found.group(1) if found else 'n/a'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
