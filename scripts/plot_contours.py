"""Plot zero-contour of phi_final.txt vs initial phi_0000.vtk (and exact shape) for an ls2d run directory.
usage: python plot_contours.py out/zalesak100 [out/vortex128 ...]
"""
import sys, os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

def read_vtk(fn):
    with open(fn) as f: lines = f.read().split("\n")
    dims = [int(x) for x in lines[4].split()[1:3]]; org = [float(x) for x in lines[5].split()[1:3]]; sp = [float(x) for x in lines[6].split()[1:3]]
    nx, ny = dims[0]-1, dims[1]-1
    i0 = lines.index("LOOKUP_TABLE default") + 1
    a = np.array([float(x) for x in lines[i0:i0+nx*ny]]).reshape(ny, nx)
    x = org[0] + (np.arange(nx)+0.5)*sp[0]; y = org[1] + (np.arange(ny)+0.5)*sp[1]
    return x, y, a

for d in sys.argv[1:]:
    x, y, p0 = read_vtk(os.path.join(d, "phi_0000.vtk"))
    x, y, p1 = read_vtk(os.path.join(d, "phi_final.vtk"))
    X, Y = np.meshgrid(x, y)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.contour(X, Y, p0, levels=[0], colors="g", linewidths=1.0)
    ax.contour(X, Y, p1, levels=[0], colors="r", linewidths=1.0)
    ax.set_aspect("equal"); ax.set_title(os.path.basename(d) + "  green=initial red=final")
    fig.savefig(os.path.join(d, "contours.png"), dpi=150, bbox_inches="tight"); print("wrote", os.path.join(d, "contours.png"))
