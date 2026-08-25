#!/usr/bin/env python3
"""Digitise the Rayleigh–Taylor interface tips from Guermond & Salgado (JCP 2009) Fig. 1 (Re = 1000).

The paper gives no table, only the density field at t_T = 1, 1.5, 1.75, 2, 2.25, 2.5 (page 11, rendered to
references/gs2009_p11-11.png at 150 dpi). Each frame is the half domain x in [0, 1/2], y in [-2, 2]; the heavy
fluid is drawn red, the light one blue. The y axis is calibrated from the coloured box (top = +2, bottom = -2,
~89 px per unit, so tip positions are accurate to about +-0.01). Spike tip = lowest red pixel, bubble tip =
highest blue pixel. Prints a table next to the ls2d 64x256 values from docs/results.md.
"""
import sys
import numpy as np
from PIL import Image

png = sys.argv[1] if len(sys.argv) > 1 else "references/gs2009_p11-11.png"
a = np.array(Image.open(png).convert("RGB")).astype(int)
R, G, B = a[..., 0], a[..., 1], a[..., 2]
red = (R >= 100) & (R > G + 60) & (R > B + 60)     # heavy fluid
blue = (B >= 100) & (B > R + 60) & (B > G + 60)    # light fluid
top, bot = 170, 560                                # Fig. 1 is the upper row of the page
col = (red[top:bot] | blue[top:bot]).sum(0) > 100
frames, s = [], None
for x in range(a.shape[1]):
    if col[x] and s is None: s = x
    if not col[x] and s is not None:
        if x - s > 20: frames.append((s, x))
        s = None
tT = [1, 1.5, 1.75, 2, 2.25, 2.5]
print("t_T   px/unit  GS spike  GS bubble")
for (x0, x1), t in zip(frames, tT):
    r, b = red[top:bot, x0:x1], blue[top:bot, x0:x1]
    rows = np.where((r | b).sum(1) > 0)[0]; yt, yb = rows.min(), rows.max()
    y = lambda row: 2 - 4 * (row - yt) / (yb - yt)
    print(f"{t:4}  {(yb - yt) / 4:6.1f}   {y(np.where(r.sum(1) > 0)[0].max()):+.3f}    {y(np.where(b.sum(1) > 0)[0].min()):+.3f}")
