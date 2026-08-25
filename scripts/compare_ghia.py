"""Compare ls2d cavity centreline profiles with Ghia, Ghia & Shin, JCP 48 (1982), Tables I and II.

usage: python3 scripts/compare_ghia.py <Re> <run dir> [<Re> <run dir> ...]
e.g.   python3 scripts/compare_ghia.py 100 out/cavity_re100_128 400 out/cavity_re400_128

Reads `references/ghia1982_cavity.txt` and, from each run directory, the `centreline_u.txt`
(u at x = 0.5) and `centreline_v.txt` (v at y = 0.5) written by ./ls2d.  The ls2d profiles are
cell-centred, so they are linearly interpolated to Ghia's tabulated locations; the two wall
locations of each table (u = 0/1 at y = 0/1, v = 0 at x = 0/1) are imposed exactly by the
boundary condition and are skipped.  Prints one markdown table per Re plus the maximum deviation.

The Re = 400 entry of Table II at x = 0.9063 (-0.23827) is a misprint in the original (see the
header of references/ghia1982_cavity.txt); it is listed but excluded from the maximum.
"""
import sys, os

BAD = {(400, "v", 0.9063)}   # printed value known to be a misprint


def read_ghia(path):
    cols = {"u": {}, "v": {}}
    which, marker = None, {"TABLE_I": "u", "TABLE_II": "v"}
    for line in open(path):
        line = line.split("#")[0].strip()
        if not line:
            continue
        if line in marker:
            which = marker[line]
            continue
        f = [float(x) for x in line.split()]
        cols[which][f[0]] = dict(zip((100, 400, 1000), f[1:]))
    return cols


def read_profile(path):
    pts = []
    for line in open(path):
        line = line.split("#")[0].strip()
        if line:
            a, b = line.split()
            pts.append((float(a), float(b)))
    return sorted(pts)


def interp(pts, q):
    if q < pts[0][0] or q > pts[-1][0]:
        return None                      # outside the cell centres: the wall value is set by the BC
    for k in range(len(pts) - 1):
        if pts[k][0] <= q <= pts[k + 1][0]:
            x0, y0 = pts[k]; x1, y1 = pts[k + 1]
            return y0 + (y1 - y0) * (q - x0) / (x1 - x0) if x1 > x0 else y0
    return pts[-1][1]


def main(argv):
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ghia = read_ghia(os.path.join(here, "references", "ghia1982_cavity.txt"))
    for re_s, d in zip(argv[0::2], argv[1::2]):
        Re = int(re_s)
        prof = {"u": read_profile(os.path.join(d, "centreline_u.txt")),
                "v": read_profile(os.path.join(d, "centreline_v.txt"))}
        worst = 0.0
        for comp, coord in (("u", "y"), ("v", "x")):
            print("\n### Re = %d, %s along the %s centreline (%s)" %
                  (Re, comp, "vertical" if comp == "u" else "horizontal", d))
            print("| %s | Ghia | ls2d | diff |\n|---|---|---|---|" % coord)
            for q in sorted(ghia[comp], reverse=True):
                ref = ghia[comp][q][Re]
                got = interp(prof[comp], q)
                if got is None:
                    print("| %.4f | %+.5f | (wall, exact) | - |" % (q, ref))
                    continue
                dev = got - ref
                flag = " *" if (Re, comp, q) in BAD else ""
                if not flag:
                    worst = max(worst, abs(dev))
                print("| %.4f | %+.5f | %+.5f | %+.5f%s |" % (q, ref, got, dev, flag))
        print("\nmax |ls2d - Ghia| at Re = %d: %.4f\n" % (Re, worst))


if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) % 2 == 0:
        sys.exit(__doc__)
    main(sys.argv[1:])
