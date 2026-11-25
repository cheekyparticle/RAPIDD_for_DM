"""Write lib/nuclear_responses/<model>.txt from the export of export_responses.wls.

    python write_responses.py <export> <model> <output> [--keep <file> --keep-A 128,129,...]

model = Anand2013 or McElvain2022. With --keep, the rows of the nuclei --keep-A are copied from <file> (a RAPIDD
table in the isospin format) instead of being taken from the export: it is how the Xe rows of McElvain2022, which
come from dmformfactor-V7.wl, are kept bit by bit (see README.md). The kept rows are compared with the exported
ones (maximum relative difference per nucleus, printed).
"""
import argparse
import math

RESPONSES = (("M", "M", "M"), ("SigmaPP", "Spp", "Spp"), ("SigmaP", "Sp", "Sp"), ("PhiPP", "Ppp", "Ppp"),
             ("PhiTildeP", "Pt", "Pt"), ("Delta", "D", "D"), ("M_PhiPP", "M", "Ppp"), ("SigmaP_Delta", "Sp", "D"))
ISOSPIN = (("0", "0"), ("1", "1"), ("0", "1"), ("1", "0"))   # ++, --, +-, -+ (order of the RAPIDD tables)
TAU = {"0": "+", "1": "-"}
ELEMENT = {6: "C", 8: "O", 9: "F", 11: "Na", 13: "Al", 14: "Si", 18: "Ar", 20: "Ca", 32: "Ge", 53: "I", 54: "Xe",
           55: "Cs"}

FORMULA = """#
# F_{O1,O2}^{tau1,tau2}(E_R) = 16 pi / (2 j_N + 1) * exp(-2y) * sum_{k=0}^{10} a_k y^k,  y = (q b / 2)^2,
#   q^2 = 2 m_T E_R, j_N the spin of the nucleus; isospin basis, tau1, tau2 = +, - (+- = (p +- n) / 2).
#   A missing row is 0. Responses and format: lib/nuclear_responses/README.md.
#"""

HEADERS = {
    "Anand2013": """# RAPIDD nuclear responses: Anand2013
# Nuclear responses of DMFormFactor v6 (N. Anand, A. L. Fitzpatrick, W. C. Haxton, arXiv:1308.6288,
#   https://www.ocf.berkeley.edu/~nanand/software/dmformfactor/) with its default one-body density matrices,
#   exported by lib/nuclear_responses/generate. DMFormFactor has no default density matrices for C, O, Al, Ar,
#   Ca and Cs. They are rounded to 5 digits, so F_M^{++}(0) differs from A^2 by up to 6e-4.""",
    "McElvain2022": """# RAPIDD nuclear responses: McElvain2022
# Nuclear responses of DMFormFactor v6 (N. Anand, A. L. Fitzpatrick, W. C. Haxton, arXiv:1308.6288) with the
#   one-body density matrices of the Elastic repository (https://github.com/Berkeley-Electroweak-Physics/Elastic;
#   W. C. Haxton, K. S. McElvain et al., arXiv:2208.07945), exported by lib/nuclear_responses/generate.
# Interaction (Elastic file) of each nucleus:
#   12C, 13C          C12_ck, C13_ck (Cohen-Kurath, p shell)
#   16O, 18O          O16_4hw, O18_usdb
#   19F               F19_usdb
#   23Na              Na23_usdb
#   27Al              Al27_usdb
#   28Si, 29Si, 30Si  Si28_usdb, Si29_usdb, Si30_usdb
#   40Ar              Ar40_Sussex
#   40Ca              Ca40_kbp
#   70-76Ge           Ge70_GCN2850 ... Ge76_GCN2850
#   127I              I127_GCN5082
#   133Cs             Cs133_GCN5082
#   128-136Xe         Xe128_GCN5082 ... Xe136_GCN5082
# usdb is the Brown-Richter 2006 refit of the USD interaction; for 16O, 4hw is used because bw is the closed core
#   without correlations. W (CaWO4) is not in the Elastic repository.
# The Xe rows come from dmformfactor-V7.wl (DMFormFactor with the Elastic GCN5082 density matrices pasted in);
#   regenerated from the Elastic files they agree to %(xe_agreement)s.""",
}


def parse_export(path):
    '''{(A, Z, response, t1, t2): (j_N, [a_0 ... a_10])} and the order of the nuclei'''
    rows, nuclei = {}, []
    for line in open(path):
        x = line.split()
        if not x:
            continue
        A, Z = int(x[0]), int(x[1])
        if (A, Z) not in nuclei:
            nuclei.append((A, Z))
        coeffs = [float(v.replace("*^", "e")) for v in x[6:]]
        assert len(coeffs) == 11, line
        rows[(A, Z, x[2], x[3], x[4])] = (float(x[5].replace("*^", "e")), coeffs)
    return rows, nuclei


def parse_table(path):
    '''RAPIDD isospin table: {(A, Z): [line tuples (A, Z, O1, O2, tau1, tau2, j_N, [a_k])]}'''
    out = {}
    for line in open(path):
        x = line.split()
        if not x or x[0].startswith("#"):
            continue
        A, Z = int(x[0]), int(x[1])
        out.setdefault((A, Z), []).append((A, Z, x[2], x[3], x[4], x[5], float(x[6]), [float(v) for v in x[7:]]))
    return out


def fmt(v):
    return "0" if v == 0 else repr(v)


def fmt_j(j):
    return "%d" % j if j == int(j) else repr(j)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("export")
    ap.add_argument("model", choices=sorted(HEADERS))
    ap.add_argument("output")
    ap.add_argument("--keep", help="RAPIDD table with the rows to keep")
    ap.add_argument("--keep-A", default="", help="comma-separated A of the nuclei whose rows are kept")
    args = ap.parse_args()

    rows, nuclei = parse_export(args.export)
    keep_A = {int(a) for a in args.keep_A.split(",") if a}
    kept = parse_table(args.keep) if args.keep else {}

    lines = []     # (A, Z, O1, O2, tau1, tau2, j_N, [a_k]) per nucleus
    agreement = {}
    for A, Z in nuclei:
        generated = []
        for name, o1, o2 in RESPONSES:
            for t1, t2 in ISOSPIN:
                j, a = rows[(A, Z, name, t1, t2)]
                if any(a):
                    generated.append((A, Z, o1, o2, TAU[t1], TAU[t2], j, a))
        if A in keep_A:
            old = kept[(A, Z)]
            gen = {r[2:6]: r[7] for r in generated}
            if {r[2:6] for r in old} != set(gen):
                raise SystemExit("kept rows of %d differ in the set of responses from the export" % A)
            worst = 0.
            for r in old:
                g = gen[r[2:6]]
                sc = max(abs(v) for v in r[7]) or 1.
                worst = max(worst, max(abs(u - v) for u, v in zip(r[7], g)) / sc)
            agreement[A] = worst
            lines.append(old)
        else:
            lines.append(generated)

    if agreement:
        print("kept rows vs export, max relative difference (scaled by the largest coefficient of the row):")
        for A, w in sorted(agreement.items()):
            print("  %d: %.1e" % (A, w))

    cells = [[(fmt(v)) for v in r[7]] for block in lines for r in block]
    width = [max(max(len(c[k]) for c in cells), len("a%d" % k)) for k in range(11)]
    jw = max(3, max(len(fmt_j(r[6])) for block in lines for r in block))
    head = HEADERS[args.model]
    if "%(xe_agreement)s" in head:
        w = max(agreement.values()) if agreement else float("nan")
        head = head % {"xe_agreement": "%.0e (relative to the largest coefficient of each row)" % w}
    out = [head, FORMULA.lstrip("\n"),
           "#" + "%4s%4s  %-5s%-5s%-6s%-6s%*s" % ("A", "Z", "O1", "O2", "tau1", "tau2", jw, "j_N")
           + "".join("  %*s" % (width[k], "a%d" % k) for k in range(11))]
    for block in lines:
        A, Z = block[0][0], block[0][1]
        out.append("")
        out.append("# %d%s (Z = %d)" % (A, ELEMENT[Z], Z))
        prev = None
        for r in block:
            if prev is not None and (r[2], r[3]) != prev:
                out.append("")
            prev = (r[2], r[3])
            out.append("%5d%4d  %-5s%-5s%-6s%-6s%*s" % (r[0], r[1], r[2], r[3], r[4], r[5], jw, fmt_j(r[6]))
                       + "".join("  %*s" % (width[k], fmt(v)) for k, v in enumerate(r[7])))
    with open(args.output, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
