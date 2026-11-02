# Python-двойник C#-решателя WellPadOptimizer (тот же алгоритм) — для проверки и расчёта результатов.
import json, math, sys, os, time
import numpy as np

DEG = math.pi / 180

class P:  # параметры (со скриншота)
    MaxOffset = 1500.0; MaxMD = 4000.0; LimitsToBottomHole = False
    Dls1 = 0.5; Dls2 = 1.5; DlsBase = 10.0
    UseAzTurn = True; MaxAzTurn = 90.0; VertLen = 40.0
    StartMaxWells = 10; EndMaxWells = 24; VariantStep = 1
    MaxHookLoad = 180.0; MaxTorque = 54.0
    CostNNS = 65.0; CostGS = 55.0
    PadCost = [(12, 16000.0), (19, 13000.0), (24, 12000.0)]
    SandCostPer1000 = 1550.0
    WellsInGroup = 4; DistInGroup = 6.0; DistBetweenGroups = 15.0
    PadWidth = 120.0; PadEndMargin = 60.0; PadMinLen = 100.0; SlopeM = 1.5
    FillHeight = {1: 1.0, 2: 2.5, 3: 4.0}
    ForbiddenBuffer = 0.0; ExtraPadGap = 0.0
    SiteStep = 100.0; SampleStep = 10.0
    Angles = [0, 30, 60, 90, 120, 150]
    # бурильная колонна для T&D
    PipeWeightAir = 0.284; MudDensity = 1.20; SteelDensity = 7.85
    FricCased = 0.25; FricOpen = 0.30; PipeRadius = 0.08
    TrajStep = 10.0
    MinSeparation = 8.0; AcMaxMD = 1000.0

def R(dls): return P.DlsBase / (dls * DEG)

def pad_build_cost(n):
    if n <= 0: return 0.0
    for lim, c in P.PadCost:
        if n <= lim: return c * n
    return math.inf

def row_len(n):
    if n <= 1: return 0.0
    g = (n - 1) // P.WellsInGroup
    return (n - 1) * P.DistInGroup + g * (P.DistBetweenGroups - P.DistInGroup)

def slot_offsets(n):
    s = [k * P.DistInGroup + (k // P.WellsInGroup) * (P.DistBetweenGroups - P.DistInGroup) for k in range(n)]
    c = (s[-1] if s else 0) / 2
    return [x - c for x in s]

def pad_len(n): return max(P.PadMinLen, row_len(n) + 2 * P.PadEndMargin)
def pad_diag(n): return math.hypot(pad_len(n), P.PadWidth)

# ---------------- карта ----------------
class Terrain:
    def __init__(self, path):
        with open(path) as f:
            hdr = {}
            for _ in range(6):
                k, v = f.readline().split(); hdr[k.lower()] = float(v)
            self.g = np.loadtxt(f, dtype=np.int16)
        self.nc, self.nr = int(hdr['ncols']), int(hdr['nrows'])
        self.x0, self.y0, self.cs = hdr['xllcorner'], hdr['yllcorner'], hdr['cellsize']
        self.nodata = int(hdr['nodata_value'])
    def cls(self, E, N):  # векторно; вне карты -> 9
        E = np.asarray(E); N = np.asarray(N)
        ix = np.floor((E - self.x0) / self.cs).astype(int)
        iy = self.nr - 1 - np.floor((N - self.y0) / self.cs).astype(int)
        ok = (ix >= 0) & (ix < self.nc) & (iy >= 0) & (iy < self.nr)
        out = np.full(E.shape, 9, dtype=np.int16)
        out[ok] = self.g[iy[ok], ix[ok]]
        out[out == self.nodata] = 9
        return out

def footprint_samples():
    """точки (s вдоль ряда, w поперёк) для макс. площадки; одинаковы для всех n"""
    half = pad_len(P.EndMaxWells) / 2 + P.ForbiddenBuffer
    hw = P.PadWidth / 2 + P.ForbiddenBuffer
    ns = int(math.ceil(half / P.SampleStep)); nw = int(math.ceil(hw / P.SampleStep))
    pts = [(i * P.SampleStep if abs(i) < ns else math.copysign(half, i), j * P.SampleStep if abs(j) < nw else math.copysign(hw, j))
           for i in range(-ns, ns + 1) for j in range(-nw, nw + 1)]
    return np.array(pts)

def sand_volume(frac, n):
    L, W = pad_len(n), P.PadWidth
    return sum(frac[c] * h * (L + P.SlopeM * h) * (W + P.SlopeM * h) for c, h in P.FillHeight.items())

def site_sand_table(terr, sites, nmax):
    """sand[site, n] (тыс.руб) и лучший угол; inf если запрещено"""
    S = footprint_samples(); ns = len(sites)
    sand = np.full((ns, nmax + 1), np.inf); best_ang = np.zeros((ns, nmax + 1), dtype=int)
    sand[:, 0] = 0
    for ai, a in enumerate(P.Angles):
        u = np.array([math.cos(a * DEG), math.sin(a * DEG)])      # ось ряда (N,E)
        v = np.array([-math.sin(a * DEG), math.cos(a * DEG)])
        dN = S[:, 0] * u[0] + S[:, 1] * v[0]; dE = S[:, 0] * u[1] + S[:, 1] * v[1]
        C = terr.cls(sites[:, 1][:, None] + dE[None, :], sites[:, 0][:, None] + dN[None, :])
        for n in range(1, nmax + 1):
            m_s = np.abs(S[:, 0]) <= pad_len(n) / 2 + P.ForbiddenBuffer + 1e-6
            m_w = np.abs(S[:, 1]) <= P.PadWidth / 2 + P.ForbiddenBuffer + 1e-6
            sub = C[:, m_s & m_w]
            forb = (sub == 9).any(axis=1)
            # доли классов — только по самой площадке (без буфера)
            m_in = (np.abs(S[:, 0]) <= pad_len(n) / 2 + 1e-6) & (np.abs(S[:, 1]) <= P.PadWidth / 2 + 1e-6)
            sub2 = C[:, m_in]; cnt = sub2.shape[1]
            vol = np.zeros(ns)
            L, W = pad_len(n), P.PadWidth
            for c, h in P.FillHeight.items():
                vol += (sub2 == c).sum(1) / cnt * h * (L + P.SlopeM * h) * (W + P.SlopeM * h)
            cost = np.where(forb, np.inf, vol / 1000 * P.SandCostPer1000)
            better = cost < sand[:, n] - 1e-9
            sand[better, n] = cost[better]; best_ang[better, n] = ai
    return sand, best_ang

# ---------------- траектория ----------------
def unit(v): return v / np.linalg.norm(v)

def solve_profile(S, T1, T3):
    """S=(N,E,0). Профиль: вертикаль -> набор1 -> тангенс -> набор2 (выход в T1 по направлению T1->T3).
       Возвращает dict или None."""
    R1, R2 = R(P.Dls1), R(P.Dls2)
    v = np.array([0, 0, 1.0]); K = np.array([S[0], S[1], P.VertLen])
    hz = np.linalg.norm(T3 - T1)
    t = unit(T3 - T1) if hz > 1e-6 else None
    d = unit(T1 - K)
    ok = False
    for it in range(300):
        th1 = math.acos(max(-1, min(1, d @ v))); P1 = K + R1 * math.tan(th1 / 2) * (v + d)
        if t is not None:
            th2 = math.acos(max(-1, min(1, d @ t)))
            if th2 > math.pi - 1e-3: return None
            P2 = T1 - R2 * math.tan(th2 / 2) * (d + t)
        else: th2 = 0.0; P2 = T1
        w = P2 - P1; Lw = np.linalg.norm(w)
        if Lw < 1e-6: return None
        dn = w / Lw
        if np.linalg.norm(dn - d) < 1e-10: ok = True; d = dn; break
        d = unit(0.5 * d + 0.5 * dn)
    if not ok: return None
    th1 = math.acos(max(-1, min(1, d @ v))); P1 = K + R1 * math.tan(th1 / 2) * (v + d)
    if t is not None:
        th2 = math.acos(max(-1, min(1, d @ t))); P2 = T1 - R2 * math.tan(th2 / 2) * (d + t)
    else: th2 = 0; P2 = T1
    L = (P2 - P1) @ d
    if L < -1e-6 or d[2] < -1e-9: return None
    md1 = P.VertLen + R1 * th1 + L + R2 * th2
    return dict(d=d, th1=th1, th2=th2, L=max(L, 0), P1=P1, P2=P2, K=K, mdT1=md1, mdT3=md1 + hz, hz=hz, t=t,
                inc_t=math.degrees(th1))

def az_of(dn, de): return (math.degrees(math.atan2(de, dn)) + 360) % 360

def az_turn(S, T1, T3):
    if np.linalg.norm(T3[:2] - T1[:2]) < 1e-6: return 0.0
    a1 = az_of(T1[0] - S[0], T1[1] - S[1]); a2 = az_of(T3[0] - T1[0], T3[1] - T1[1])
    d = abs(a1 - a2) % 360
    return min(d, 360 - d)

def pair_check(S, T1, T3):
    """полная проверка ограничений, возвращает (prof|None, причина)"""
    off = math.hypot((T3 if P.LimitsToBottomHole else T1)[0] - S[0], (T3 if P.LimitsToBottomHole else T1)[1] - S[1])
    if off > P.MaxOffset + 1e-6: return None, 'offset'
    if P.UseAzTurn and az_turn(S, T1, T3) > P.MaxAzTurn + 1e-6: return None, 'azturn'
    pr = solve_profile(S, T1, T3)
    if pr is None: return None, 'profile'
    md = pr['mdT3'] if P.LimitsToBottomHole else pr['mdT1']
    if md > P.MaxMD + 1e-6: return None, 'md'
    return pr, ''

def drill_cost(pr): return P.CostNNS * pr['mdT1'] + P.CostGS * pr['hz']

def arc_points(P0, a, b, Rr, md0, step):
    th = math.acos(max(-1, min(1, a @ b)))
    out = []
    if th < 1e-9: return out, P0, md0
    n = unit(b - (a @ b) * a)
    Larc = Rr * th; k = 1
    while k * step < Larc - 1e-6:
        ph = k * step / Rr
        out.append((md0 + k * step, P0 + Rr * (math.sin(ph) * a + (1 - math.cos(ph)) * n), math.cos(ph) * a + math.sin(ph) * n))
        k += 1
    Pend = P0 + Rr * (math.sin(th) * a + (1 - math.cos(th)) * n)
    out.append((md0 + Larc, Pend, b))
    return out, Pend, md0 + Larc

def line_points(P0, d, L, md0, step):
    out = []; k = 1
    while k * step < L - 1e-6:
        out.append((md0 + k * step, P0 + k * step * d, d)); k += 1
    if L > 1e-9: out.append((md0 + L, P0 + L * d, d))
    return out, P0 + L * d, md0 + L

def survey(S, T1, T3, pr, step=None):
    step = step or P.TrajStep
    v = np.array([0, 0, 1.0]); S3 = np.array([S[0], S[1], 0.0])
    pts = [(0.0, S3, v)]
    seg, Pc, md = line_points(S3, v, P.VertLen, 0.0, step); pts += seg
    seg, Pc, md = arc_points(Pc, v, pr['d'], R(P.Dls1), md, step); pts += seg
    seg, Pc, md = line_points(Pc, pr['d'], pr['L'], md, step); pts += seg
    if pr['t'] is not None:
        seg, Pc, md = arc_points(Pc, pr['d'], pr['t'], R(P.Dls2), md, step); pts += seg
        seg, Pc, md = line_points(Pc, pr['t'], pr['hz'], md, step); pts += seg
    rows = []
    prev = None
    for md, p, tg in pts:
        inc = math.degrees(math.acos(max(-1, min(1, tg[2]))))
        az = az_of(tg[0], tg[1]) if inc > 1e-6 else 0.0
        dls = 0.0
        if prev is not None and md - prev[0] > 1e-6:
            ang = math.degrees(math.acos(max(-1, min(1, prev[2] @ tg))))
            dls = ang / (md - prev[0]) * P.DlsBase
        rows.append([md, inc, az, p[0], p[1], p[2], dls])
        prev = (md, p, tg)
    return np.array(rows)

def torque_drag(sv, mdT1):
    """мягкая модель Johancsik: нагрузка на крюке при подъёме (тс) и момент при вращении над забоем (кН·м)"""
    bf = 1 - P.MudDensity / P.SteelDensity; w = P.PipeWeightAir * bf
    Tp = 0.0; Tr = 0.0; Tq = 0.0
    for i in range(len(sv) - 1, 0, -1):
        md2, inc2, az2 = sv[i, 0], sv[i, 1] * DEG, sv[i, 2] * DEG
        md1, inc1, az1 = sv[i - 1, 0], sv[i - 1, 1] * DEG, sv[i - 1, 2] * DEG
        dL = md2 - md1
        if dL <= 0: continue
        ia = (inc1 + inc2) / 2; di = inc2 - inc1
        da = (az2 - az1 + math.pi) % (2 * math.pi) - math.pi
        if sv[i - 1, 1] < 1e-6 or sv[i, 1] < 1e-6: da = 0.0
        mu = P.FricOpen if md1 >= mdT1 else P.FricCased
        Np = math.hypot(Tp * da * math.sin(ia), Tp * di + w * dL * math.sin(ia))
        Tp = Tp + w * dL * math.cos(ia) + mu * Np
        Nr = math.hypot(Tr * da * math.sin(ia), Tr * di + w * dL * math.sin(ia))
        Tr = Tr + w * dL * math.cos(ia)
        Tq += mu * Nr * P.PipeRadius
    return Tp / 9.80665, Tq

# ---------------- оптимизация ----------------
class Problem:
    def __init__(self, inp, terr):
        self.raw = json.load(open(inp)); self.terr = terr
        r = self.raw['massiveDataVariants'][0]['trajBuildRequest']
        self.wells = []
        for i, w in enumerate(r):
            T1 = np.array([w['T1'][0], w['T1'][1], w['T1'][2]]); T3 = np.array([w['T3'][0], w['T3'][1], w['T3'][2]])
            self.wells.append(dict(i=i, uid=w['uid'], name=f"W{i+1:02d}", T1=T1, T3=T3))
        self.nw = len(self.wells)
        self.real = self.raw.get('trajRealWells', [])

    def build_sites(self):
        T1 = np.array([w['T1'][:2] for w in self.wells])
        lo = T1.min(0) - P.MaxOffset; hi = T1.max(0) + P.MaxOffset
        lo = np.floor(lo / P.SiteStep) * P.SiteStep; hi = np.ceil(hi / P.SiteStep) * P.SiteStep
        nn = np.arange(lo[0], hi[0] + 1, P.SiteStep); ee = np.arange(lo[1], hi[1] + 1, P.SiteStep)
        G = np.array([(n, e) for n in nn for e in ee])
        dmin = np.min(np.hypot(G[:, None, 0] - T1[None, :, 0], G[:, None, 1] - T1[None, :, 1]), axis=1)
        self.sites = G[dmin <= P.MaxOffset]
        print('sites', len(self.sites))

    def build_matrix(self):
        ns = len(self.sites); C = np.full((ns, self.nw), np.inf)
        reasons = {}
        for j, w in enumerate(self.wells):
            T1, T3 = w['T1'], w['T3']
            off = np.hypot(self.sites[:, 0] - T1[0], self.sites[:, 1] - T1[1])
            for s in np.nonzero(off <= P.MaxOffset)[0]:
                S = (self.sites[s, 0], self.sites[s, 1])
                pr, why = pair_check(S, T1, T3)
                if pr is not None: C[s, j] = drill_cost(pr)
                else: reasons[why] = reasons.get(why, 0) + 1
        self.C = C
        print('feasible pairs', np.isfinite(C).sum(), reasons)

def total_cost(sol, prob, sand):
    tot = 0.0
    for p in sol:
        n = len(p['w']); tot += pad_build_cost(n) + sand[p['s'], n] + prob.C[p['s'], p['w']].sum()
    return tot

def pad_fixed(s, n, sand): return pad_build_cost(n) + sand[s, n] if n > 0 else 0.0

def optimize(prob, sand, K, max_pads=None):
    C = prob.C; ns, nw = C.shape; sites = prob.sites
    Dmin = pad_diag(K) + P.ExtraPadGap
    def site_ok(s, pads, skip=None):
        for q, p in enumerate(pads):
            if q == skip: continue
            if math.hypot(sites[s, 0] - sites[p['s'], 0], sites[s, 1] - sites[p['s'], 1]) < Dmin: return False
        return True
    # 1) жадная инициализация по мин. средней стоимости на скважину
    pads = []; un = np.ones(nw, bool)
    while un.any():
        Cu = np.where(un[None, :], C, np.inf)
        order = np.argsort(Cu, axis=1); Cs = np.take_along_axis(Cu, order, 1)
        best = (math.inf, -1, 0)
        cum = np.cumsum(Cs[:, :K], axis=1)
        for k in range(1, min(K, un.sum()) + 1):
            avg = (np.array([pad_build_cost(k)]) + sand[:, k] + cum[:, k - 1]) / k
            for s in np.argsort(avg)[:50]:
                if not np.isfinite(avg[s]): break
                if avg[s] >= best[0]: break
                if site_ok(s, pads): best = (avg[s], s, k); break
        if best[1] < 0:
            raise RuntimeError('не удалось разместить скважины: ' + str(np.nonzero(un)[0]))
        s, k = best[1], best[2]
        ws = list(order[s, :k]); pads.append(dict(s=s, w=ws)); un[ws] = False
    def relocate(pads, pi):
        p = pads[pi]; n = len(p['w'])
        cost = sand[:, n] + C[:, p['w']].sum(1)
        curc = cost[p['s']]
        for s in np.argsort(cost, kind='stable')[:300]:
            if cost[s] >= curc - 1e-6: return False
            if site_ok(s, pads, pi): p['s'] = s; return True
        return False

    def close_pad(pads, pi):
        """закрыть площадку pi, скважины раздать остальным; None если невозможно"""
        trial = [dict(s=p['s'], w=list(p['w'])) for p in pads]
        rem = trial.pop(pi)
        for w in sorted(rem['w'], key=lambda w: (np.min(C[[t['s'] for t in trial], w]), w)):
            cand = [(C[t['s'], w] + pad_fixed(t['s'], len(t['w']) + 1, sand) - pad_fixed(t['s'], len(t['w']), sand), ti)
                    for ti, t in enumerate(trial) if len(t['w']) < K and np.isfinite(C[t['s'], w])]
            if not cand: return None
            trial[min(cand)[1]]['w'].append(w)
        for ti in range(len(trial)): relocate(trial, ti)
        return trial

    def local_search(pads):
        cur = total_cost(pads, prob, sand)
        for it in range(200):
            # a) перенос скважины в другой куст
            for pi in range(len(pads)):
                for w in list(pads[pi]['w']):
                    p = pads[pi]; n_p = len(p['w'])
                    if n_p == 0: break
                    base_p = pad_fixed(p['s'], n_p, sand) - pad_fixed(p['s'], n_p - 1, sand)
                    best = (-1e-6, -1)
                    for qi, q in enumerate(pads):
                        if qi == pi or not q['w'] or len(q['w']) >= K or not np.isfinite(C[q['s'], w]): continue
                        n_q = len(q['w'])
                        delta = pad_fixed(q['s'], n_q + 1, sand) - pad_fixed(q['s'], n_q, sand) - base_p + C[q['s'], w] - C[p['s'], w]
                        if delta < best[0]: best = (delta, qi)
                    if best[1] >= 0:
                        p['w'].remove(w); pads[best[1]]['w'].append(w)
            pads = [p for p in pads if p['w']]
            # b) обмен скважинами между кустами
            for pi in range(len(pads)):
                for qi in range(pi + 1, len(pads)):
                    p, q = pads[pi], pads[qi]
                    for a in list(p['w']):
                        for b in list(q['w']):
                            if a not in p['w'] or b not in q['w']: continue
                            dl = C[q['s'], a] + C[p['s'], b] - C[p['s'], a] - C[q['s'], b]
                            if np.isfinite(dl) and dl < -1e-6:
                                p['w'].remove(a); q['w'].remove(b); p['w'].append(b); q['w'].append(a)
            # c) перенос площадок
            for pi in range(len(pads)): relocate(pads, pi)
            # d) закрытие куста
            if len(pads) > 1:
                for pi in range(len(pads)):
                    trial = close_pad(pads, pi)
                    if trial is not None and total_cost(trial, prob, sand) < total_cost(pads, prob, sand) - 1e-6:
                        pads = trial; break
            newc = total_cost(pads, prob, sand)
            if newc > cur - 1e-6: cur = newc; break
            cur = newc
        return pads, cur

    pads, cur = local_search(pads)
    # 3) принудительное сокращение числа кустов (для сравнения вариантов)
    while max_pads is not None and len(pads) > max_pads:
        best = None
        for pi in range(len(pads)):
            trial = close_pad(pads, pi)
            if trial is None: continue
            c = total_cost(trial, prob, sand)
            if best is None or c < best[0]: best = (c, trial)
        if best is None: return None, math.inf
        pads, cur = local_search(best[1])
        pads = pads if len(pads) <= len(best[1]) else pads
    return pads, cur

# ---------------- финальная оценка варианта ----------------
def pad_geometry(terr, center, n, mean_hz_az):
    """подбор ориентации для фактической площадки"""
    S = footprint_samples(); best = None
    for a in P.Angles:
        u = np.array([math.cos(a * DEG), math.sin(a * DEG)]); v = np.array([-math.sin(a * DEG), math.cos(a * DEG)])
        m_b = (np.abs(S[:, 0]) <= pad_len(n) / 2 + P.ForbiddenBuffer + 1e-6) & (np.abs(S[:, 1]) <= P.PadWidth / 2 + P.ForbiddenBuffer + 1e-6)
        m_in = (np.abs(S[:, 0]) <= pad_len(n) / 2 + 1e-6) & (np.abs(S[:, 1]) <= P.PadWidth / 2 + 1e-6)
        dN = S[:, 0] * u[0] + S[:, 1] * v[0]; dE = S[:, 0] * u[1] + S[:, 1] * v[1]
        cl = terr.cls(center[1] + dE, center[0] + dN)
        forb = bool((cl[m_b] == 9).any())
        sub = cl[m_in]; frac = {c: float((sub == c).mean()) for c in P.FillHeight}
        vol = sand_volume(frac, n); cost = vol / 1000 * P.SandCostPer1000
        perp = abs(((a - (mean_hz_az + 90)) + 90) % 180 - 90) if mean_hz_az is not None else 0
        key = (forb, round(cost, 0), perp)
        if best is None or key < best[0]: best = (key, a, vol, cost, frac, forb)
    return best

def evaluate(prob, pads_in, name, kind, K=None, fixed_wh=None):
    """pads_in: list of dict(center=(N,E), wells=[idx]); fixed_wh: dict well->(N,E) для эталонных вариантов"""
    res = dict(name=name, kind=kind, K=K, pads=[], wells=[])
    for pi, p in enumerate(pads_in):
        ws = p['wells']; n = len(ws); c = np.array(p['center'])
        hz_az = [az_of(*(prob.wells[w]['T3'][:2] - prob.wells[w]['T1'][:2])) for w in ws if prob.wells[w]['T3'][:2].tolist() != prob.wells[w]['T1'][:2].tolist()]
        mean_az = None
        if hz_az:
            sx = sum(math.cos(2 * a * DEG) for a in hz_az); sy = sum(math.sin(2 * a * DEG) for a in hz_az)
            mean_az = (math.degrees(math.atan2(sy, sx)) / 2) % 180
        if fixed_wh is not None:
            # ориентация по устьям (PCA)
            pts = np.array([fixed_wh[w] for w in ws]); dd = pts - pts.mean(0)
            if n > 1:
                ev = np.linalg.eigh(dd.T @ dd)[1][:, -1]; ang = az_of(ev[0], ev[1]) % 180
            else: ang = 0
            _, _, vol, scost, frac, forb = pad_geometry_fixed(prob.terr, c, n, ang)
            angle = ang; slots = {w: fixed_wh[w] for w in ws}
        else:
            (forb, _, _), angle, vol, scost, frac, forb = pad_geometry(prob.terr, c, n, mean_az)
            u = np.array([math.cos(angle * DEG), math.sin(angle * DEG)])
            offs = slot_offsets(n)
            order = sorted(ws, key=lambda w: (prob.wells[w]['T1'][:2] - c) @ u)
            order = improve_slots(prob, order, c, u, offs)
            slots = {w: tuple(c + offs[k] * u) for k, w in enumerate(order)}
        pad = dict(id=pi + 1, name=f"К-{pi+1}", N=float(c[0]), E=float(c[1]), n=n, angle=float(angle),
                   length=pad_len(n), width=P.PadWidth, terrain={str(k): round(v, 3) for k, v in frac.items()},
                   sand_volume=vol, sand_cost=scost, build_cost=pad_build_cost(n), forbidden=forb, wells=[])
        for k, w in enumerate(list(slots)):
            W = prob.wells[w]; S = slots[w]
            pr = solve_profile(S, W['T1'], W['T3'])
            off = math.hypot(W['T1'][0] - S[0], W['T1'][1] - S[1])
            turn = az_turn(S, W['T1'], W['T3'])
            viol = []
            if off > P.MaxOffset + 1e-6: viol.append('отход')
            if P.UseAzTurn and turn > P.MaxAzTurn + 1e-6: viol.append('разворот')
            rec = dict(well=W['name'], uid=W['uid'], pad=pad['name'], whN=float(S[0]), whE=float(S[1]), offset=off, azturn=turn,
                       T1=W['T1'].tolist(), T3=W['T3'].tolist(), hz=float(np.linalg.norm(W['T3'] - W['T1'])))
            if pr is None:
                viol.append('профиль'); rec.update(mdT1=None, mdT3=None, cost=None, hook=None, torque=None, maxinc=None, traj=[])
            else:
                sv = survey(S, W['T1'], W['T3'], pr)
                hook, tq = torque_drag(sv, pr['mdT1'])
                md = pr['mdT3'] if P.LimitsToBottomHole else pr['mdT1']
                if md > P.MaxMD + 1e-6: viol.append('MD')
                if hook > P.MaxHookLoad: viol.append('крюк')
                if tq > P.MaxTorque: viol.append('момент')
                rec.update(mdT1=pr['mdT1'], mdT3=pr['mdT3'], tangent=pr['L'], inc_tangent=pr['inc_t'],
                           az_tangent=az_of(pr['d'][0], pr['d'][1]), cost=drill_cost(pr), cost_nns=P.CostNNS * pr['mdT1'],
                           cost_gs=P.CostGS * pr['hz'], hook=hook, torque=tq, maxdls=float(sv[:, 6].max()), survey=sv)
            rec['viol'] = viol
            pad['wells'].append(rec['well']); res['wells'].append(rec)
        res['pads'].append(pad)
    # антиколлизия внутри куста
    by_pad = {}
    for r in res['wells']: by_pad.setdefault(r['pad'], []).append(r)
    for r in res['wells']: r['min_sep'] = None
    for lst in by_pad.values():
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                a, b = lst[i].get('survey'), lst[j].get('survey')
                if a is None or b is None: continue
                A = a[a[:, 0] > P.VertLen][:, 3:6]; B = b[b[:, 0] > P.VertLen][:, 3:6]
                if len(A) == 0 or len(B) == 0: continue
                dm = float(np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)).min())
                for r in (lst[i], lst[j]):
                    r['min_sep'] = dm if r['min_sep'] is None else min(r['min_sep'], dm)
    # итоги
    wl = res['wells']
    res['n_pads'] = len(res['pads']); res['n_wells'] = len(wl)
    res['md_total'] = sum(r['mdT3'] or 0 for r in wl)
    res['drill_cost'] = sum(r['cost'] or 0 for r in wl)
    res['pad_build_cost'] = sum(p['build_cost'] for p in res['pads'])
    res['sand_cost'] = sum(p['sand_cost'] for p in res['pads'])
    res['sand_volume'] = sum(p['sand_volume'] for p in res['pads'])
    res['total_cost'] = res['drill_cost'] + res['pad_build_cost'] + res['sand_cost']
    res['violations'] = sum(1 for r in wl if r['viol']) + sum(1 for p in res['pads'] if p['forbidden'])
    return res

def shallow_path(prob, w, S):
    W = prob.wells[w]; pr = solve_profile(S, W['T1'], W['T3'])
    if pr is None: return None, math.inf
    sv = survey(S, W['T1'], W['T3'], pr)
    sv = sv[(sv[:, 0] > P.VertLen) & (sv[:, 0] <= P.AcMaxMD)]
    return sv[:, 3:6], drill_cost(pr)

def min_sep(a, b):
    if a is None or b is None or len(a) == 0 or len(b) == 0: return math.inf
    m = min(len(a), len(b))  # сравнение по одинаковым MD (одинаковая сетка MD)
    return float(np.sqrt(((a[:m] - b[:m]) ** 2).sum(1)).min())

def improve_slots(prob, order, c, u, offs):
    """перестановки скважин по слотам: максимизация мин. сближения (до AcMaxMD), затем мин. стоимость"""
    n = len(order)
    if n < 2: return order
    cache = {}
    def path(w, k):
        if (w, k) not in cache: cache[(w, k)] = shallow_path(prob, w, tuple(c + offs[k] * u))
        return cache[(w, k)]
    def score(ordr):
        ms = math.inf; cost = 0.0
        for k in range(n):
            pa, cc = path(ordr[k], k); cost += cc
            for j in range(k + 1, n):
                ms = min(ms, min_sep(pa, path(ordr[j], j)[0]))
        return (min(ms, P.MinSeparation), -cost)
    cur = score(order); order = list(order)
    for _ in range(4):
        imp = False
        for i in range(n):
            for j in range(i + 1, n):
                order[i], order[j] = order[j], order[i]
                sc = score(order)
                if sc > cur: cur = sc; imp = True
                else: order[i], order[j] = order[j], order[i]
        if not imp: break
    return order

def pad_geometry_fixed(terr, center, n, ang):
    S = footprint_samples(); a = ang
    u = np.array([math.cos(a * DEG), math.sin(a * DEG)]); v = np.array([-math.sin(a * DEG), math.cos(a * DEG)])
    m_in = (np.abs(S[:, 0]) <= pad_len(n) / 2 + 1e-6) & (np.abs(S[:, 1]) <= P.PadWidth / 2 + 1e-6)
    dN = S[:, 0] * u[0] + S[:, 1] * v[0]; dE = S[:, 0] * u[1] + S[:, 1] * v[1]
    cl = terr.cls(center[1] + dE, center[0] + dN)
    forb = bool((cl[m_in] == 9).any())
    sub = cl[m_in]; frac = {c: float((sub == c).mean()) for c in P.FillHeight}
    vol = sand_volume(frac, n)
    return None, a, vol, vol / 1000 * P.SandCostPer1000, frac, forb

def main(inp, mapf, outdir):
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    terr = Terrain(mapf); prob = Problem(inp, terr)
    prob.build_sites(); prob.build_matrix()
    print('matrix', round(time.time() - t0, 1), 's')
    sand, _ = site_sand_table(terr, prob.sites, P.EndMaxWells)
    print('sand', round(time.time() - t0, 1), 's', np.isfinite(sand[:, 1]).mean())
    results = []
    # эталонные варианты из входного файла
    for vi, v in enumerate(prob.raw['massiveDataVariants']):
        r = v['trajBuildRequest']; uid2i = {w['uid']: w['i'] for w in prob.wells}
        key = {(round(w['T1'][0], 2), round(w['T1'][1], 2), round(w['T3'][0], 2), round(w['T3'][1], 2)): w['i'] for w in prob.wells}
        groups = {}; wh = {}
        for w in r:
            i = key[(round(w['T1'][0], 2), round(w['T1'][1], 2), round(w['T3'][0], 2), round(w['T3'][1], 2))]
            groups.setdefault(w['pad_id'], []).append(i); wh[i] = (w['wh_north'], w['wh_east'])
        pads_in = [dict(center=tuple(np.mean([wh[i] for i in g], 0)), wells=g) for _, g in sorted(groups.items())]
        res = evaluate(prob, pads_in, f"Исходный вариант {vi}", 'input', fixed_wh=wh)
        results.append(res); print(res['name'], round(res['total_cost'] / 1000, 1), 'млн', res['violations'])
    for sc_name, use_az in [('Разворот ≤%g°' % P.MaxAzTurn, True), ('Без ограничения разворота', False)]:
        P.UseAzTurn = use_az
        if not use_az: prob.build_matrix()
        prev = None
        for K in range(P.StartMaxWells, P.EndMaxWells + 1, P.VariantStep):
            pads, cost = optimize(prob, sand, K)
            sig = sorted((int(p['s']), tuple(sorted(p['w']))) for p in pads)
            if prev is not None and sig == prev[0]:
                res = dict(prev[1]); res = {k: v for k, v in prev[1].items()}
                res['name'] = f"{sc_name}; ≤{K} скв./куст"; res['K'] = K; res['same_as'] = prev[1]['name']
                results.append(res); print(res['name'], '= same'); continue
            pads_in = [dict(center=tuple(prob.sites[p['s']]), wells=sorted(p['w'])) for p in pads]
            res = evaluate(prob, pads_in, f"{sc_name}; ≤{K} скв./куст", 'opt', K=K)
            res['opt_estimate'] = cost; res['scenario'] = sc_name; res['same_as'] = None
            results.append(res); prev = (sig, res)
            print(res['name'], res['n_pads'], [p['n'] for p in res['pads']], round(cost / 1000, 1), round(res['total_cost'] / 1000, 1), 'млн', res['violations'], round(time.time() - t0, 1))
        P.UseAzTurn = True
    # варианты с ограничением числа кустов (K = макс.), сценарий с разворотом
    P.UseAzTurn = False; prob.build_matrix()
    base = [r for r in results if r.get('scenario', '').startswith('Без') and not r.get('same_as')][-1]
    for mp in range(base['n_pads'] - 1, 0, -1):
        pads, cost = optimize(prob, sand, P.EndMaxWells, max_pads=mp)
        if pads is None: print('max_pads', mp, 'infeasible'); break
        P.UseAzTurn = False
        pads_in = [dict(center=tuple(prob.sites[p['s']]), wells=sorted(p['w'])) for p in pads]
        res = evaluate(prob, pads_in, f"Без ограничения разворота; не более {mp} кустов", 'opt', K=P.EndMaxWells)
        res['opt_estimate'] = cost; res['scenario'] = 'Лимит кустов'; res['same_as'] = None
        results.append(res)
        print(res['name'], res['n_pads'], [p['n'] for p in res['pads']], round(cost / 1000, 1), round(res['total_cost'] / 1000, 1), 'млн', res['violations'], round(time.time() - t0, 1))
    P.UseAzTurn = True
    return prob, results

if __name__ == '__main__':
    prob, results = main(sys.argv[1], sys.argv[2], sys.argv[3])
    import pickle; pickle.dump((results,), open(os.path.join(sys.argv[3], 'results.pkl'), 'wb'))
