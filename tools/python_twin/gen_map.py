# Генератор синтетической карты желательности (ESRI ASCII Grid).
# Коды: 1 - суходол, 2 - болото, 3 - водоём, 9 - запрещено.
import json, sys, numpy as np
inp, out = sys.argv[1], sys.argv[2]
d = json.load(open(inp))
r = d['massiveDataVariants'][0]['trajBuildRequest']
E = [w['T1'][1] for w in r] + [w['T3'][1] for w in r]
N = [w['T1'][0] for w in r] + [w['T3'][0] for w in r]
cell = 50.0; m = 2500
x0 = np.floor((min(E) - m) / 1000) * 1000; y0 = np.floor((min(N) - m) / 1000) * 1000
x1 = np.ceil((max(E) + m) / 1000) * 1000;  y1 = np.ceil((max(N) + m) / 1000) * 1000
nc = int((x1 - x0) / cell); nr = int((y1 - y0) / cell)
rng = np.random.default_rng(20261102)
def octave(scale):
    gy, gx = int(nr * cell / scale) + 3, int(nc * cell / scale) + 3
    g = rng.random((gy, gx))
    yy = (np.arange(nr) + 0.5) * cell / scale; xx = (np.arange(nc) + 0.5) * cell / scale
    iy = yy.astype(int); ix = xx.astype(int); fy = (yy - iy)[:, None]; fx = (xx - ix)[None, :]
    fy = fy * fy * (3 - 2 * fy); fx = fx * fx * (3 - 2 * fx)
    a = g[iy][:, ix]; b = g[iy][:, ix + 1]; c = g[iy + 1][:, ix]; dd = g[iy + 1][:, ix + 1]
    return (a * (1 - fx) + b * fx) * (1 - fy) + (c * (1 - fx) + dd * fx) * fy
f = octave(3000) * 0.5 + octave(1200) * 0.3 + octave(500) * 0.15 + octave(200) * 0.05
q55, q88 = np.quantile(f, [0.50, 0.88])
cls = np.where(f > q88, 3, np.where(f > q55, 2, 1)).astype(int)
# Запретные зоны: ООПТ/инфраструктура (круги) и коридор ЛЭП/трубопровода.
yc = y1 - (np.arange(nr) + 0.5) * cell   # строки сверху вниз
xc = x0 + (np.arange(nc) + 0.5) * cell
X, Y = np.meshgrid(xc, yc)
t1 = np.array([[w['T1'][1], w['T1'][0]] for w in r])
for k in range(14):
    c = t1[rng.integers(len(t1))] + rng.normal(0, 900, 2)
    rad = rng.uniform(250, 550)
    cls[(X - c[0]) ** 2 + (Y - c[1]) ** 2 < rad ** 2] = 9
# два линейных коридора шириной 120 м
for (ax, ay, bx, by) in [(x0, 6771200, x1, 6773900), (x0, 6813800, x1, 6811300)]:
    dx, dy = bx - ax, by - ay; L = np.hypot(dx, dy)
    dist = np.abs((X - ax) * dy - (Y - ay) * dx) / L
    cls[dist < 60] = 9
with open(out, 'w') as fo:
    fo.write(f"ncols {nc}\nnrows {nr}\nxllcorner {x0:.1f}\nyllcorner {y0:.1f}\ncellsize {cell}\nNODATA_value -9999\n")
    for row in cls: fo.write(' '.join(map(str, row)) + '\n')
u, c = np.unique(cls, return_counts=True); print(nc, nr, x0, y0, dict(zip(u.tolist(), (c / c.sum()).round(3).tolist())))
