namespace WellPadOptimizer;

/// <summary>Куст в ходе оптимизации: индекс площадки-кандидата (−1 — не размещён) и список скважин.</summary>
public sealed class PadAssignment
{
    public int Site = -1;
    public List<int> Wells = new();
    public PadAssignment Clone() => new() { Site = Site, Wells = new List<int>(Wells) };
}

/// <summary>
/// Размещение кустов и распределение скважин.
/// Целевая функция: Σ кусты [строительство(n) + отсыпка(место, n)] + Σ скважины [бурение(место, скважина) + штраф].
///  1. Кандидатные площадки — сетка SiteGridStep в радиусе SoftRadiusFactor·MaxOffset от T1.
///  2. Матрица C[место, скважина] = бурение + штраф за каждое нарушенное ограничение (отход, разворот, MD).
///     Штраф на порядки больше стоимости бурения, поэтому оптимизатор сначала минимизирует число нарушений,
///     затем их величину и только потом затраты. Если профиль не строится — ∞.
///  3. Таблица отсыпки Sand[место, n] — по карте желательности, лучшая ориентация, ∞ если задета запретная зона.
///  4. Жёсткое ограничение размера куста: MinWellsPerPad ≤ n ≤ K.
///     Скважины делятся на группы (участки); для каждой группы перебирается число кустов p
///     (ceil(m/K) ≤ p ≤ floor(m/Nmin)), строится разбиение k-средних с исправлением размеров,
///     затем локальный поиск: перенос скважины, обмен, перенос площадки. Лучшие p по группам комбинируются.
/// </summary>
public sealed class Optimizer
{
    private readonly Parameters _p;
    private readonly InputData _data;
    private readonly TerrainMap _map;
    public double[] SiteN, SiteE;
    public double[,] C;           // тыс. руб. (с учётом штрафов)
    public double[,] Sand;        // тыс. руб.
    public int NSites => SiteN.Length;
    public int NWells => _data.Wells.Count;
    public int FeasiblePairs { get; private set; }
    public int CleanPairs { get; private set; }

    public Optimizer(Parameters p, InputData data, TerrainMap map) { _p = p; _data = data; _map = map; }

    private double SoftRadius => _p.MaxOffset * _p.SoftRadiusFactor;

    // -------------------------------------------------------------- подготовка

    public void BuildSites()
    {
        var w = _data.Wells;
        double step = _p.SiteGridStep, r = SoftRadius;
        double n0 = Math.Floor((w.Min(x => x.T1.N) - r) / step) * step;
        double n1 = Math.Ceiling((w.Max(x => x.T1.N) + r) / step) * step;
        double e0 = Math.Floor((w.Min(x => x.T1.E) - r) / step) * step;
        double e1 = Math.Ceiling((w.Max(x => x.T1.E) + r) / step) * step;
        var ln = new List<double>(); var le = new List<double>();
        for (double n = n0; n <= n1 + 1e-6; n += step)
            for (double e = e0; e <= e1 + 1e-6; e += step)
            {
                bool near = false;
                foreach (var t in w)
                    if (Geo.Hypot(t.T1.N - n, t.T1.E - e) <= r) { near = true; break; }
                if (near) { ln.Add(n); le.Add(e); }
            }
        SiteN = ln.ToArray(); SiteE = le.ToArray();
    }

    /// <summary>Стоимость пары с штрафами; nViol — число нарушенных ограничений.</summary>
    public static double SoftCost(Parameters p, Vec3 wh, WellTarget w, out int nViol)
    {
        nViol = 0;
        var prof = TrajectoryCalculator.Solve(p, wh, w.T1, w.T3);
        if (prof == null) { nViol = 9; return double.PositiveInfinity; }
        double pen = 0;
        var lim = p.LimitsToBottomHole ? w.T3 : w.T1;
        double ex = Geo.Hypot(lim.N - wh.N, lim.E - wh.E) - p.MaxOffset;
        if (ex > 1e-6) { pen += p.PenaltyViolation + p.PenaltyOffsetPerM * ex; nViol++; }
        if (p.UseAzimuthTurn)
        {
            ex = TrajectoryCalculator.AzimuthTurn(wh, w.T1, w.T3) - p.MaxAzimuthTurn;
            if (ex > 1e-6) { pen += p.PenaltyViolation + p.PenaltyAzimuthPerDeg * ex; nViol++; }
        }
        ex = (p.LimitsToBottomHole ? prof.MdT3 : prof.MdT1) - p.MaxMD;
        if (ex > 1e-6) { pen += p.PenaltyViolation + p.PenaltyMdPerM * ex; nViol++; }
        return TrajectoryCalculator.DrillCost(p, prof) + pen;
    }

    public void BuildCostMatrix()
    {
        int ns = NSites, nw = NWells;
        C = new double[ns, nw];
        int feas = 0, clean = 0;
        double r = SoftRadius;
        Parallel.For(0, nw, j =>
        {
            var w = _data.Wells[j];
            int f = 0, c = 0;
            for (int s = 0; s < ns; s++)
            {
                C[s, j] = double.PositiveInfinity;
                if (Geo.Hypot(SiteN[s] - w.T1.N, SiteE[s] - w.T1.E) > r) continue;
                C[s, j] = SoftCost(_p, new Vec3(SiteN[s], SiteE[s], 0), w, out int nv);
                if (!double.IsInfinity(C[s, j])) { f++; if (nv == 0) c++; }
            }
            Interlocked.Add(ref feas, f); Interlocked.Add(ref clean, c);
        });
        FeasiblePairs = feas; CleanPairs = clean;
    }

    /// <summary>Точки опробования площадки (s — вдоль ряда устьев, w — поперёк) для максимального куста.</summary>
    public static List<(double s, double w)> FootprintSamples(Parameters p, int nmax)
    {
        double half = p.PadLength(nmax) / 2 + p.ForbiddenBuffer, hw = p.PadWidth / 2 + p.ForbiddenBuffer;
        int ns = (int)Math.Ceiling(half / p.FootprintSampleStep), nw = (int)Math.Ceiling(hw / p.FootprintSampleStep);
        var pts = new List<(double, double)>();
        for (int i = -ns; i <= ns; i++)
            for (int j = -nw; j <= nw; j++)
                pts.Add((Math.Abs(i) < ns ? i * p.FootprintSampleStep : Math.CopySign(half, i),
                         Math.Abs(j) < nw ? j * p.FootprintSampleStep : Math.CopySign(hw, j)));
        return pts;
    }

    /// <summary>Оценка площадки в точке: запрещена ли, доли классов, объём и стоимость отсыпки.</summary>
    public static (bool forbidden, Dictionary<int, double> frac, double volume, double cost) EvaluateFootprint(
        Parameters p, TerrainMap map, List<(double s, double w)> samples, double cn, double ce, int n, double angleDeg)
    {
        double a = angleDeg * Geo.Deg;
        double uN = Math.Cos(a), uE = Math.Sin(a), vN = -Math.Sin(a), vE = Math.Cos(a);
        double halfL = p.PadLength(n) / 2, halfW = p.PadWidth / 2;
        bool forb = false; int cnt = 0;
        var counts = new Dictionary<int, int>();
        foreach (var (s, w) in samples)
        {
            bool inBuf = Math.Abs(s) <= halfL + p.ForbiddenBuffer + 1e-6 && Math.Abs(w) <= halfW + p.ForbiddenBuffer + 1e-6;
            if (!inBuf) continue;
            int cls = map.ClassAt(ce + s * uE + w * vE, cn + s * uN + w * vN);
            if (cls == p.ForbiddenCode) forb = true;
            if (Math.Abs(s) <= halfL + 1e-6 && Math.Abs(w) <= halfW + 1e-6)
            {
                cnt++;
                counts[cls] = counts.TryGetValue(cls, out int c) ? c + 1 : 1;
            }
        }
        var frac = new Dictionary<int, double>();
        foreach (var t in p.Terrain) frac[t.Code] = cnt > 0 && counts.TryGetValue(t.Code, out int c) ? (double)c / cnt : 0;
        double vol = p.SandVolume(frac, n);
        return (forb, frac, vol, p.SandCost(vol));
    }

    public void BuildSandTable()
    {
        int nmax = _p.MaxPadCapacity;
        Sand = new double[NSites, nmax + 1];
        var samples = FootprintSamples(_p, nmax);
        Parallel.For(0, NSites, s =>
        {
            Sand[s, 0] = 0;
            for (int n = 1; n <= nmax; n++)
            {
                double best = double.PositiveInfinity;
                foreach (var ang in _p.PadOrientations)
                {
                    var r = EvaluateFootprint(_p, _map, samples, SiteN[s], SiteE[s], n, ang);
                    if (!r.forbidden && r.cost < best) best = r.cost;
                }
                Sand[s, n] = best;
            }
        });
    }

    // -------------------------------------------------------------- целевая функция

    private double PadFixed(int s, int n) => n <= 0 ? 0 : _p.PadBuildCost(n) + Sand[s, n];

    public double TotalCost(List<PadAssignment> pads)
    {
        double t = 0;
        foreach (var p in pads)
        {
            if (p.Site < 0) return double.PositiveInfinity;
            t += PadFixed(p.Site, p.Wells.Count);
            foreach (var w in p.Wells) t += C[p.Site, w];
        }
        return t;
    }

    // -------------------------------------------------------------- группы скважин

    /// <summary>Группы скважин, которые могут обслуживаться общими площадками (T1 ближе 2·R).</summary>
    public List<List<int>> Groups()
    {
        var w = _data.Wells; int n = w.Count; double lim = 2 * SoftRadius;
        var lab = Enumerable.Range(0, n).ToArray();
        int Find(int x) { while (lab[x] != x) x = lab[x] = lab[lab[x]]; return x; }
        for (int i = 0; i < n; i++)
            for (int k = i + 1; k < n; k++)
                if (Geo.Hypot(w[i].T1.N - w[k].T1.N, w[i].T1.E - w[k].T1.E) < lim) lab[Find(i)] = Find(k);
        return Enumerable.Range(0, n).GroupBy(Find).Select(g => g.ToList()).OrderBy(g => g.Min()).ToList();
    }

    /// <summary>Диапазон общего числа кустов при размере куста [Nmin, K].</summary>
    public (int min, int max) PadCountRange(int K)
    {
        int lo = 0, hi = 0;
        foreach (var g in Groups()) { lo += (g.Count + K - 1) / K; hi += g.Count / _p.MinWellsPerPad; }
        return (lo, hi);
    }

    // -------------------------------------------------------------- оптимизация

    /// <summary>
    /// Оптимизация при размере куста [MinWellsPerPad, K]. nPads — заданное общее число кустов (или null).
    /// Возвращает (null, ∞), если скважины нельзя разбить на кусты такого размера.
    /// </summary>
    public (List<PadAssignment> pads, double cost) Optimize(int K, int? nPads = null)
    {
        var table = new List<Dictionary<int, (List<PadAssignment> pads, double cost)>>();
        foreach (var g in Groups())
        {
            var opts = new Dictionary<int, (List<PadAssignment>, double)>();
            int pMin = (g.Count + K - 1) / K, pMax = g.Count / _p.MinWellsPerPad;
            for (int p = pMin; p <= pMax; p++)
            {
                var (pads, c) = OptimizeGroup(g, K, p);
                if (pads != null && !double.IsInfinity(c)) opts[p] = (pads, c);
            }
            if (opts.Count == 0) return (null, double.PositiveInfinity);
            table.Add(opts);
        }
        // перебор комбинаций числа кустов по группам
        var combos = new List<(int[] ps, double cost, int n)> { (Array.Empty<int>(), 0.0, 0) };
        foreach (var opts in table)
        {
            var next = new List<(int[], double, int)>();
            foreach (var cb in combos)
                foreach (var kv in opts)
                    next.Add((cb.ps.Append(kv.Key).ToArray(), cb.cost + kv.Value.cost, cb.n + kv.Key));
            combos = next;
        }
        if (nPads.HasValue) combos = combos.Where(x => x.n == nPads.Value).ToList();
        if (combos.Count == 0) return (null, double.PositiveInfinity);
        var best = combos.OrderBy(x => x.cost).First();
        var result = new List<PadAssignment>();
        for (int i = 0; i < table.Count; i++) result.AddRange(table[i][best.ps[i]].pads);
        return (result, best.cost);
    }

    private (List<PadAssignment>, double) OptimizeGroup(List<int> wells, int K, int p)
    {
        int m = wells.Count, nmin = _p.MinWellsPerPad, ns = NSites;
        if (p * nmin > m || p * K < m) return (null, double.PositiveInfinity);
        double dmin = _p.PadDiagonal(K) + _p.ExtraPadGap;
        var tn = wells.Select(w => _data.Wells[w].T1.N).ToArray();
        var te = wells.Select(w => _data.Wells[w].T1.E).ToArray();

        bool SiteOk(int s, List<PadAssignment> others, int skip)
        {
            for (int q = 0; q < others.Count; q++)
            {
                if (q == skip || others[q].Site < 0) continue;
                int o = others[q].Site;
                if (Geo.Hypot(SiteN[s] - SiteN[o], SiteE[s] - SiteE[o]) < dmin) return false;
            }
            return true;
        }

        // лучшая площадка для набора скважин куста (только улучшение, если площадка уже есть)
        bool Place(List<PadAssignment> list, int pi)
        {
            var pd = list[pi]; int n = pd.Wells.Count;
            var cost = new double[ns]; var idx = new int[ns];
            for (int s = 0; s < ns; s++)
            {
                double c = Sand[s, n];
                foreach (var w in pd.Wells) c += C[s, w];
                cost[s] = c; idx[s] = s;
            }
            Array.Sort((double[])cost.Clone(), idx);
            for (int i = 0; i < Math.Min(500, ns); i++)
            {
                int s = idx[i];
                if (double.IsInfinity(cost[s])) break;
                if (pd.Site >= 0 && cost[s] >= cost[pd.Site] - 1e-6) return false;
                if (SiteOk(s, list, pi)) { pd.Site = s; return true; }
            }
            return false;
        }

        // k-средних по T1 с исправлением размеров кластеров до [nmin, K]
        List<PadAssignment> KMeans(int seed)
        {
            var cn = new List<double> { tn[seed] }; var ce = new List<double> { te[seed] };
            while (cn.Count < p)
            {
                int far = 0; double fd = -1;
                for (int i = 0; i < m; i++)
                {
                    double d = double.PositiveInfinity;
                    for (int k = 0; k < cn.Count; k++) d = Math.Min(d, Geo.Hypot(tn[i] - cn[k], te[i] - ce[k]));
                    if (d > fd) { fd = d; far = i; }
                }
                cn.Add(tn[far]); ce.Add(te[far]);
            }
            var lab = new int[m];
            for (int it = 0; it < 30; it++)
            {
                for (int i = 0; i < m; i++)
                {
                    int bk = 0; double bd = double.PositiveInfinity;
                    for (int k = 0; k < p; k++) { double d = Geo.Hypot(tn[i] - cn[k], te[i] - ce[k]); if (d < bd) { bd = d; bk = k; } }
                    lab[i] = bk;
                }
                bool moved = false;
                for (int k = 0; k < p; k++)
                {
                    var mem = Enumerable.Range(0, m).Where(i => lab[i] == k).ToList();
                    if (mem.Count == 0) continue;
                    double nn = mem.Average(i => tn[i]), ee = mem.Average(i => te[i]);
                    if (Math.Abs(nn - cn[k]) > 1e-6 || Math.Abs(ee - ce[k]) > 1e-6) moved = true;
                    cn[k] = nn; ce[k] = ee;
                }
                if (!moved) break;
            }
            double D(int i, int k) => Geo.Hypot(tn[i] - cn[k], te[i] - ce[k]);
            for (int guard = 0; guard < 10 * m; guard++)
            {
                var sz = new int[p]; foreach (var l in lab) sz[l]++;
                int big = Array.FindIndex(sz, x => x > K), small = Array.FindIndex(sz, x => x < nmin);
                if (big < 0 && small < 0) break;
                double bestD = double.PositiveInfinity; int bi = -1, bq = -1;
                if (big >= 0)
                {
                    for (int i = 0; i < m; i++)
                        if (lab[i] == big)
                            for (int q = 0; q < p; q++)
                                if (q != big && sz[q] < K && D(i, q) - D(i, big) < bestD) { bestD = D(i, q) - D(i, big); bi = i; bq = q; }
                }
                else
                {
                    for (int i = 0; i < m; i++)
                        if (lab[i] != small && sz[lab[i]] > nmin && D(i, small) - D(i, lab[i]) < bestD) { bestD = D(i, small) - D(i, lab[i]); bi = i; bq = small; }
                }
                if (bi < 0) return null;
                lab[bi] = bq;
            }
            var res = new List<PadAssignment>();
            for (int k = 0; k < p; k++) res.Add(new PadAssignment { Wells = Enumerable.Range(0, m).Where(i => lab[i] == k).Select(i => wells[i]).ToList() });
            return res;
        }

        (List<PadAssignment>, double) LocalSearch(List<PadAssignment> list)
        {
            double cur = TotalCost(list);
            for (int it = 0; it < 100; it++)
            {
                // a) перенос скважины в другой куст (с сохранением размеров [nmin, K])
                for (int pi = 0; pi < list.Count; pi++)
                {
                    foreach (var w in list[pi].Wells.ToList())
                    {
                        var pd = list[pi]; int np = pd.Wells.Count;
                        if (np - 1 < nmin) break;
                        double baseP = PadFixed(pd.Site, np) - PadFixed(pd.Site, np - 1);
                        double best = -1e-6; int bq = -1;
                        for (int qi = 0; qi < list.Count; qi++)
                        {
                            var q = list[qi];
                            if (qi == pi || q.Wells.Count >= K || double.IsInfinity(C[q.Site, w])) continue;
                            int nq = q.Wells.Count;
                            double delta = PadFixed(q.Site, nq + 1) - PadFixed(q.Site, nq) - baseP + C[q.Site, w] - C[pd.Site, w];
                            if (delta < best) { best = delta; bq = qi; }
                        }
                        if (bq >= 0) { pd.Wells.Remove(w); list[bq].Wells.Add(w); }
                    }
                }
                // b) обмен скважинами
                for (int pi = 0; pi < list.Count; pi++)
                    for (int qi = pi + 1; qi < list.Count; qi++)
                    {
                        var a = list[pi]; var b = list[qi];
                        foreach (var x in a.Wells.ToList())
                            foreach (var y in b.Wells.ToList())
                            {
                                if (!a.Wells.Contains(x) || !b.Wells.Contains(y)) continue;
                                double dl = C[b.Site, x] + C[a.Site, y] - C[a.Site, x] - C[b.Site, y];
                                if (!double.IsNaN(dl) && !double.IsInfinity(dl) && dl < -1e-6)
                                {
                                    a.Wells.Remove(x); b.Wells.Remove(y); a.Wells.Add(y); b.Wells.Add(x);
                                }
                            }
                    }
                // c) перенос площадок
                for (int pi = 0; pi < list.Count; pi++) Place(list, pi);
                double nc = TotalCost(list);
                if (nc > cur - 1e-6) { cur = nc; break; }
                cur = nc;
            }
            return (list, cur);
        }

        // старты: скважины, наиболее удалённые от центра группы
        double mn = tn.Average(), me = te.Average();
        var starts = p > 1
            ? Enumerable.Range(0, m).OrderByDescending(i => Geo.Hypot(tn[i] - mn, te[i] - me)).Take(6).ToList()
            : new List<int> { 0 };
        List<PadAssignment> bestPads = null; double bestCost = double.PositiveInfinity;
        foreach (var st in starts)
        {
            var pads = KMeans(st);
            if (pads == null) continue;
            bool ok = true;
            for (int pi = 0; pi < pads.Count && ok; pi++) ok = Place(pads, pi);
            if (!ok) continue;
            var (sol, c) = LocalSearch(pads);
            if (c < bestCost - 1e-6) { bestCost = c; bestPads = sol; }
        }
        return (bestPads, bestCost);
    }
}
