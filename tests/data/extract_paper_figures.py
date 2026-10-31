"""Extract the data of the vector figures of arXiv:1912.08204 (Catena, Emken, Spaldin,
Tarantino) used by tests/test_electron_paper.py, and write them as CSV files in
tests/data/paper_1912.08204/.

The figures are Mathematica PDFs: the curves are polylines whose vertices are the plotted
points, the bars of Total_Rate_Comparison are filled rectangles. The axes are calibrated
from the tick marks and their labels (log scales, checked to be linear in log10 within
0.002 decades). Points on or outside the plot frame (clipped) are dropped.

Usage (requires PyMuPDF):
    python tests/data/extract_paper_figures.py [path/to/arXiv-1912.08204v2/figures] [--plot DIR]
--plot DIR writes a PNG per figure with the extracted data, to be compared by eye with
the original figure.
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "paper_1912.08204")
DEFAULT_FIGURES = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", "..", "..", "..",
                                                "arXiv-1912.08204v2", "figures"))

# colors of the paper (rounded to 2 digits)
BLACK = (0.0, 0.0, 0.0)
BLUE = (0.25, 0.46, 0.89)
X10 = (0.48, 0.77, 0.98)
X1T = (0.24, 0.33, 0.83)
DS50 = (0.77, 0.5, 0.81)
DS50_B = (0.69, 0.33, 0.75)
SOLID, DASHED, DOTTED, DASHDOT = "[] 0", "[ 4 4 ] 0", "[ 1 4 ] 0", "[ 1 4 4 4 ] 0"

SPEC_UNITS = ("E_e [eV]", "dR_ion/dlnE_e [kg^-1 yr^-1]")
LIM_C = ("m_chi [MeV]", "c_i^s")
LIM_L = ("m_chi [MeV]", "c_i^l")

# figure -> (units, {(color, dashes): curve name}); spectra: Xenon, 5 outer shells
FIGURES = {
    "Xe_coup_13_C": (SPEC_UNITS, {(BLUE, SOLID): "m1000MeV", (BLUE, DASHED): "m100MeV", (BLUE, DOTTED): "m10MeV"}),
    "Xe_coup_13_L": (SPEC_UNITS, {(BLUE, SOLID): "m1000MeV", (BLUE, DASHED): "m100MeV", (BLUE, DOTTED): "m10MeV"}),
    "Xe_coup_7_C_mDM_10": (SPEC_UNITS, {(BLACK, SOLID): "total", (BLUE, SOLID): "res1",
                                        (BLUE, DASHED): "res2", (BLUE, DOTTED): "res3"}),
    # res. 4 is negative (dagger): its absolute value is plotted
    "Xe_coup_15_L_mDM_100": (SPEC_UNITS, {(BLACK, SOLID): "total", (BLUE, SOLID): "res1",
                                          (BLUE, DASHED): "res3", (BLUE, DOTTED): "res4_negative"}),
    "Spectrum_Anapole": (SPEC_UNITS, {(BLACK, SOLID): "total", (BLUE, SOLID): "O8", (BLUE, DASHED): "O9"}),
    "Spectrum_Magnetic_Dipole": (SPEC_UNITS, {(BLACK, SOLID): "total", (BLUE, SOLID): "O1", (BLUE, DASHED): "O4",
                                              (BLUE, DOTTED): "O5", (BLUE, DASHDOT): "O6"}),
    "Spectrum_Electric_Dipole": (SPEC_UNITS, {(BLACK, SOLID): "total", (BLUE, SOLID): "O11"}),
    "Constraint_Coupling_C": (LIM_C, {(c, d): "%s_c%s" % (e, o)
                                      for c, e in ((X10, "XENON10"), (X1T, "XENON1T"), (DS50, "DarkSide50"))
                                      for d, o in ((SOLID, "1"), (DASHED, "8"), (DOTTED, "13"))}),
    "Constraint_Coupling_L": (LIM_L, {(c, d): "%s_c%s" % (e, o)
                                      for c, e in ((X10, "XENON10"), (X1T, "XENON1T"), (DS50, "DarkSide50"))
                                      for d, o in ((SOLID, "1"), (DASHED, "8"), (DOTTED, "13"))}),
    "Constraint_Anapole": (("m_chi [MeV]", "g/Lambda^2 [GeV^-2]"),
                           {(X10, SOLID): "XENON10", (X1T, SOLID): "XENON1T", (DS50_B, SOLID): "DarkSide50"}),
    "Constraint_Magnetic_Dipole": (("m_chi [MeV]", "g/Lambda [GeV^-1]"),
                                   {(X10, SOLID): "XENON10", (X1T, SOLID): "XENON1T", (DS50_B, SOLID): "DarkSide50"}),
    "Constraint_Electric_Dipole": (("m_chi [MeV]", "g/Lambda [GeV^-1]"),
                                   {(X10, SOLID): "XENON10", (X1T, SOLID): "XENON1T", (DS50_B, SOLID): "DarkSide50"}),
}
TOTAL_RATE = "Total_Rate_Comparison"
TOTAL_RATE_OPERATORS = (1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
TOTAL_RATE_COLORS = {(0.62, 0.2, 0.7): "Ar", BLUE: "Xe"}


def rcolor(c):
    return tuple(round(v, 2) for v in c) if c else None


def segments(drawings):
    '''Thin black straight lines (axes frame and ticks) as (x0, y0, x1, y1)'''
    segs = []
    for d in drawings:
        if d["type"] != "s" or rcolor(d.get("color")) != BLACK or (d.get("width") or 0) > 0.5:
            continue
        for it in d["items"]:
            if it[0] == "l":
                segs.append((it[1].x, it[1].y, it[2].x, it[2].y))
    return segs


def find_frames(segs, min_len=100.):
    '''Plot frames (x0, y0, x1, y1) from the long horizontal and vertical thin black lines'''
    hor = sorted({(round(min(a, c), 1), round(max(a, c), 1), round(b, 1)) for a, b, c, d in segs
                  if abs(b - d) < 0.01 and abs(a - c) >= min_len})
    ver = sorted({(round(a, 1), round(min(b, d), 1), round(max(b, d), 1)) for a, b, c, d in segs
                  if abs(a - c) < 0.01 and abs(b - d) >= min_len})
    frames = []
    for x0, x1, yb in hor:
        for x0b, x1b, yt in hor:
            if yt < yb and abs(x0 - x0b) < 0.2 and abs(x1 - x1b) < 0.2:
                if any(abs(vx - x0) < 0.2 and abs(vy0 - yt) < 0.2 and abs(vy1 - yb) < 0.2 for vx, vy0, vy1 in ver):
                    frames.append((x0, yt, x1, yb))
    # keep the minimal frames (no other frame inside)
    frames = sorted(set(frames))
    return [f for f in frames if not any(g != f and g[0] >= f[0] - 0.2 and g[2] <= f[2] + 0.2
                                         and g[1] >= f[1] - 0.2 and g[3] <= f[3] + 0.2 for g in frames)]


def ticks(segs, frame, axis):
    '''Positions of the major (longest) tick marks on the bottom (axis "x") or left (axis "y")
    side of frame'''
    x0, y0, x1, y1 = frame
    pos = {}
    for a, b, c, d in segs:
        if axis == "x" and abs(a - c) < 0.01 and x0 - 0.2 <= a <= x1 + 0.2 and abs(max(b, d) - y1) < 0.2 \
                and 0 < abs(b - d) < 15:
            pos[round(a, 2)] = abs(b - d)
        if axis == "y" and abs(b - d) < 0.01 and y0 - 0.2 <= b <= y1 + 0.2 and abs(min(a, c) - x0) < 0.2 \
                and 0 < abs(a - c) < 15:
            pos[round(b, 2)] = abs(a - c)
    if not pos:
        return np.array([])
    lmax = max(pos.values())
    return np.array(sorted(x for x, l in pos.items() if l > 0.8 * lmax))


def text_spans(page):
    out = []
    for b in page.get_text("dict")["blocks"]:
        for line in b.get("lines", []):
            for s in line["spans"]:
                if s["text"].strip():
                    out.append((s["text"].strip(), pymupdf.Rect(s["bbox"]), s["size"]))
    return out


def tick_labels(spans, frame, axis):
    '''Tick labels as (position, value): "10" + superscript exponent, or a plain number'''
    x0, y0, x1, y1 = frame
    if axis == "x":
        cand = [s for s in spans if y1 < s[1].y0 < y1 + 20 and x0 - 25 < s[1].x0 < x1 + 25]
    else:
        cand = [s for s in spans if x0 - 20 < s[1].x1 < x0 and y0 - 8 < s[1].y0 and s[1].y1 < y1 + 8]
    if not cand:
        return []
    size = max(s[2] for s in cand)
    groups = []
    for b in (s for s in cand if s[2] > 0.9 * size):
        # superscript exponent: smaller spans chained to the right of the base, raised
        g = [b]
        for s in sorted((s for s in cand if s[2] < 0.9 * size), key=lambda s: s[1].x0):
            if abs(s[1].x0 - g[-1][1].x1) < 1.5 and b[1].y0 - 6 < s[1].y0 and s[1].y1 < b[1].y1:
                g.append(s)
        groups.append(g)
    labels = []
    for g in groups:
        base, exp = g[0][0], "".join(s[0] for s in g[1:]).replace("−", "-")
        if not re.fullmatch(r"\d+(\.\d+)?", base) or (exp and not re.fullmatch(r"-?\d+", exp)):
            continue
        value = float(base) ** int(exp) if exp else float(base)
        if exp and base != "10":
            continue
        r = g[0][1]
        for s in g[1:]:
            r = r | s[1]
        pos = (r.x0 + r.x1) / 2 if axis == "x" else (g[0][1].y0 + g[0][1].y1) / 2
        labels.append((pos, value))
    return labels


def calibrate(segs, spans, frame, axis):
    '''log10(value) = a + b * coordinate from the labelled ticks'''
    tk = ticks(segs, frame, axis)
    labels = tick_labels(spans, frame, axis)
    pts = []
    for pos, value in labels:
        if len(tk) and np.min(np.abs(tk - pos)) < 4.:
            pts.append((tk[np.argmin(np.abs(tk - pos))], np.log10(value)))
    pts = np.array(sorted(set(pts)))
    if len(pts) < 2:
        raise RuntimeError("cannot calibrate the %s axis of frame %s" % (axis, frame))
    b, a = np.polyfit(pts[:, 0], pts[:, 1], 1)
    res = np.max(np.abs(a + b * pts[:, 0] - pts[:, 1]))
    if res > 2e-3:
        raise RuntimeError("the %s axis is not logarithmic (residual %.3g decades)" % (axis, res))
    return a, b


def path_points(d):
    pts = []
    for it in d["items"]:
        if it[0] == "l":
            seq = (it[1], it[2])
        elif it[0] == "c":
            seq = (it[1], it[4])
        else:
            continue
        for p in seq:
            if not pts or abs(pts[-1][0] - p.x) > 1e-6 or abs(pts[-1][1] - p.y) > 1e-6:
                pts.append((p.x, p.y))
    return np.array(pts)


def extract_curves(page, name, styles):
    drawings = page.get_drawings()
    segs = segments(drawings)
    spans = text_spans(page)
    frames = find_frames(segs)
    if len(frames) != 1:
        raise RuntimeError("%s: expected one plot frame, found %s" % (name, frames))
    frame = frames[0]
    ax, bx = calibrate(segs, spans, frame, "x")
    ay, by = calibrate(segs, spans, frame, "y")
    curves = {}
    for d in drawings:
        if d["type"] != "s" or len(d["items"]) < 8:
            continue
        key = (rcolor(d.get("color")), d.get("dashes"))
        if key not in styles:
            continue
        if styles[key] in curves:
            raise RuntimeError("%s: more than one curve with style %s" % (name, key))
        p = path_points(d)
        x0, y0, x1, y1 = frame
        inside = (p[:, 0] > x0 + 0.05) & (p[:, 0] < x1 - 0.05) & (p[:, 1] > y0 + 0.05) & (p[:, 1] < y1 - 0.05)
        p = p[inside]
        curves[styles[key]] = np.column_stack((10**(ax + bx * p[:, 0]), 10**(ay + by * p[:, 1])))
    missing = set(styles.values()) - set(curves)
    if missing:
        raise RuntimeError("%s: curves not found: %s" % (name, missing))
    return curves


def extract_total_rate(page):
    '''Bars of Total_Rate_Comparison: rows (interaction, operator, target, mDM [MeV], R/c^2 [kg^-1 yr^-1])'''
    drawings = page.get_drawings()
    segs = segments(drawings)
    spans = text_spans(page)
    frames = sorted(find_frames(segs), key=lambda f: f[1])
    if len(frames) != 2:
        raise RuntimeError("Total_Rate_Comparison: expected two panels, found %s" % frames)
    rows = []
    for frame, interaction in zip(frames, ("contact", "long-range")):
        ay, by = calibrate(segs, spans, frame, "y")
        bars = defaultdict(list)
        for d in drawings:
            target = TOTAL_RATE_COLORS.get(rcolor(d.get("fill")))
            if d["type"] not in ("f", "fs") or target is None:
                continue
            r = d["rect"]
            if frame[1] - 1 <= r.y0 and r.y1 <= frame[3] + 1:
                bars[target].append(((r.x0 + r.x1) / 2, r.y0))
        for target, bb in bars.items():
            bb.sort()
            if len(bb) != 3 * len(TOTAL_RATE_OPERATORS):
                raise RuntimeError("Total_Rate_Comparison: %d bars for %s" % (len(bb), target))
            for i, (xc, ytop) in enumerate(bb):
                rows.append((interaction, TOTAL_RATE_OPERATORS[i // 3], target, (10, 100, 1000)[i % 3],
                             10**(ay + by * ytop)))
    return rows


def write_csv(path, header, cols):
    with open(path, "w", newline="") as f:
        f.write("# %s\n" % header)
        w = csv.writer(f)
        for row in zip(*cols):
            w.writerow(["%.6g" % v if isinstance(v, float) else v for v in row])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("figures", nargs="?", default=DEFAULT_FIGURES)
    parser.add_argument("--plot", default=None)
    args = parser.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    extracted = {}
    for name, (units, styles) in FIGURES.items():
        page = pymupdf.open(os.path.join(args.figures, name + ".pdf"))[0]
        curves = extract_curves(page, name, styles)
        extracted[name] = (units, curves)
        for curve, xy in curves.items():
            write_csv(os.path.join(OUT_DIR, "%s__%s.csv" % (name, curve)),
                      "arXiv:1912.08204 figures/%s.pdf, curve %s; columns: %s, %s" % (name, curve, *units),
                      (xy[:, 0], xy[:, 1]))
            print("%-28s %-18s %3d points" % (name, curve, len(xy)))
    page = pymupdf.open(os.path.join(args.figures, TOTAL_RATE + ".pdf"))[0]
    rows = extract_total_rate(page)
    write_csv(os.path.join(OUT_DIR, TOTAL_RATE + ".csv"),
              "arXiv:1912.08204 figures/%s.pdf; columns: interaction, operator, target, m_chi [MeV], "
              "R_ion/c_i^2 [kg^-1 yr^-1]" % TOTAL_RATE, list(zip(*rows)))
    print("%-28s %d bars" % (TOTAL_RATE, len(rows)))

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        os.makedirs(args.plot, exist_ok=True)
        for name, (units, curves) in extracted.items():
            fig, ax = plt.subplots(figsize=(5, 3.5))
            for curve, xy in curves.items():
                ax.loglog(xy[:, 0], xy[:, 1], ".-", ms=3, label=curve)
            ax.set_xlabel(units[0])
            ax.set_ylabel(units[1])
            ax.set_title(name, fontsize=9)
            ax.legend(fontsize=6)
            fig.tight_layout()
            fig.savefig(os.path.join(args.plot, name + "_extracted.png"), dpi=110)
            plt.close(fig)
        fig, axes = plt.subplots(2, 1, figsize=(10, 6))
        for ax, interaction in zip(axes, ("contact", "long-range")):
            for target, color in (("Ar", "purple"), ("Xe", "royalblue")):
                sel = [r for r in rows if r[0] == interaction and r[2] == target]
                xs = [TOTAL_RATE_OPERATORS.index(r[1]) * 4 + (0, 1, 2)[(10, 100, 1000).index(r[3])]
                      + (0.4 if target == "Xe" else 0) for r in sel]
                ax.bar(xs, [r[4] for r in sel], width=0.4, color=color, log=True)
            ax.set_xticks([4 * i + 1 for i in range(len(TOTAL_RATE_OPERATORS))])
            ax.set_xticklabels(["O%d" % o for o in TOTAL_RATE_OPERATORS])
        fig.savefig(os.path.join(args.plot, TOTAL_RATE + "_extracted.png"), dpi=110)
    return 0


if __name__ == "__main__":
    sys.exit(main())
