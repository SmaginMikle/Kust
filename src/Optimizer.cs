namespace WellPadOptimizer;

/// <summary>Куст в ходе оптимизации: индекс площадки-кандидата и список скважин.</summary>
public sealed class PadAssignment
{
    public int Site;
    public List<int> Wells = new();
    public PadAssignment Clone() => new() { Site = Site, Wells = new List<int>(Wells) };
}

/// <summary>
/// Размещение кустов и распределение скважин.
/// Целевая функция: Σ кусты [строительство(n) + отсыпка(место, n)] + Σ скважины [бурение(место, скважина)].
///  1. Кандидатные площадки — регулярная сетка (SiteGridStep) в радиусе MaxOffset от T1.
///  2. Матрица C[место, скважина] — стоимость бурения при выполнении всех ограничений (иначе ∞).
///  3. Таблица отсыпки Sand[место, n] — по карте желательности, лучшая ориентация, ∞ если задета запретная зона.
///  4. Жадная инициализация (минимум средней стоимости на скважину) + локальный поиск:
///     перенос скважины, обмен скважинами, перенос площадки, закрытие куста.
/// </summary>
public sealed class Optimizer
{
    private readonly Parameters _p;
    private readonly InputData _data;
    private readonly TerrainMap _map;
    public double[] SiteN, SiteE;
    public double[,] C;           // тыс. руб.
    public double[,] Sand;        // тыс. руб.
    public int NSites => SiteN.Length;
    public int NWells => _data.Wells.Count;
    public Dictionary<CheckResult, int> RejectStats { get; } = new();

    public Optimizer(Parameters p, InputData data, TerrainMap map) { _p = p; _data = data; _map = map; }

    // -------------------------------------------------------------- подготовка

    public void BuildSites()
    {
        var w = _data.Wells;
        double step = _p.SiteGridStep;
        double n0 = Math.Floor((w.Min(x => x.T1.N) - _p.MaxOffset) / step) * step;
        double n1 = Math.Ceiling((w.Max(x => x.T1.N) + _p.MaxOffset) / step) * step;
        double e0 = Math.Floor((w.Min(x => x.T1.E) - _p.MaxOffset) / step) * step;
        double e1 = Math.Ceiling((w.Max(x => x.T1.E) + _p.MaxOffset) / step) * step;
        var ln = new List<double>(); var le = new List<double>();
        for (double n = n0; n <= n1 + 1e-6; n += step)
            for (double e = e0; e <= e1 + 1e-6; e += step)
            {
                bool near = false;
                foreach (var t in w)
                    if (Geo.Hypot(t.T1.N - n, t.T1.E - e) <= _p.MaxOffset) { near = true; break; }
                if (near) { ln.Add(n); le.Add(e); }
            }
        SiteN = ln.ToArray(); SiteE = le.ToArray();
    }

    public void BuildCostMatrix()
    {
        int ns = NSites, nw = NWells;
        C = new double[ns, nw];
        var stats = new int[Enum.GetValues<CheckResult>().Length];
        Parallel.For(0, nw, j =>
        {
            var w = _data.Wells[j];
            var local = new int[stats.Length];
            for (int s = 0; s < ns; s++)
            {
                C[s, j] = double.PositiveInfinity;
                if (Geo.Hypot(SiteN[s] - w.T1.N, SiteE[s] - w.T1.E) > _p.MaxOffset) continue;
                var (prof, res) = TrajectoryCalculator.Check(_p, new Vec3(SiteN[s], SiteE[s], 0), w);
                local[(int)res]++;
                if (prof != null) C[s, j] = TrajectoryCalculator.DrillCost(_p, prof);
            }
            lock (stats) for (int k = 0; k < stats.Length; k++) stats[k] += local[k];
        });
        RejectStats.Clear();
        foreach (CheckResult r in Enum.GetValues<CheckResult>()) RejectStats[r] = stats[(int)r];
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
            t += PadFixed(p.Site, p.Wells.Count);
            foreach (var w in p.Wells) t += C[p.Site, w];
        }
        return t;
    }

    // -------------------------------------------------------------- оптимизация

    /// <summary>Оптимизация при максимальном числе скважин на кусте K и (опционально) лимите числа кустов.</summary>
    public (List<PadAssignment> pads, double cost) Optimize(int K, int? maxPads = null)
    {
        double dmin = _p.PadDiagonal(K) + _p.ExtraPadGap;
        int nw = NWells, ns = NSites;

        bool SiteOk(int s, List<PadAssignment> others, int skip = -1)
        {
            for (int q = 0; q < others.Count; q++)
            {
                if (q == skip) continue;
                int o = others[q].Site;
                if (Geo.Hypot(SiteN[s] - SiteN[o], SiteE[s] - SiteE[o]) < dmin) return false;
            }
            return true;
        }

        // 1) жадная инициализация
        var pads = new List<PadAssignment>();
        var unassigned = new bool[nw]; Array.Fill(unassigned, true);
        int left = nw;
        while (left > 0)
        {
            var cand = new List<(double avg, int s, int k)>();
            var bestPerSite = new int[ns][];
            for (int s = 0; s < ns; s++)
            {
                var lst = new List<(double c, int w)>();
                for (int w = 0; w < nw; w++)
                    if (unassigned[w] && !double.IsInfinity(C[s, w])) lst.Add((C[s, w], w));
                if (lst.Count == 0) continue;
                lst.Sort((x, y) => x.c != y.c ? x.c.CompareTo(y.c) : x.w.CompareTo(y.w));
                bestPerSite[s] = lst.Select(x => x.w).ToArray();
                double cum = 0;
                for (int k = 1; k <= Math.Min(K, lst.Count); k++)
                {
                    cum += lst[k - 1].c;
                    double avg = (PadFixed(s, k) + cum) / k;
                    if (!double.IsInfinity(avg)) cand.Add((avg, s, k));
                }
            }
            cand.Sort((x, y) => x.avg != y.avg ? x.avg.CompareTo(y.avg) : x.s != y.s ? x.s.CompareTo(y.s) : x.k.CompareTo(y.k));
            int chosen = -1;
            for (int i = 0; i < cand.Count; i++)
                if (SiteOk(cand[i].s, pads)) { chosen = i; break; }
            if (chosen < 0)
                throw new InvalidOperationException("Не удалось разместить скважины: " +
                    string.Join(", ", Enumerable.Range(0, nw).Where(w => unassigned[w]).Select(w => _data.Wells[w].Name)));
            var (_, cs, ck) = cand[chosen];
            var pa = new PadAssignment { Site = cs };
            for (int i = 0; i < ck; i++) { int w = bestPerSite[cs][i]; pa.Wells.Add(w); unassigned[w] = false; left--; }
            pads.Add(pa);
        }

        bool Relocate(List<PadAssignment> list, int pi)
        {
            var p = list[pi]; int n = p.Wells.Count;
            var cost = new double[ns]; var idx = new int[ns];
            for (int s = 0; s < ns; s++)
            {
                double c = Sand[s, n];
                foreach (var w in p.Wells) c += C[s, w];
                cost[s] = c; idx[s] = s;
            }
            double cur = cost[p.Site];
            Array.Sort((double[])cost.Clone(), idx);
            int lim = Math.Min(300, ns);
            for (int i = 0; i < lim; i++)
            {
                int s = idx[i];
                if (cost[s] >= cur - 1e-6) return false;
                if (SiteOk(s, list, pi)) { p.Site = s; return true; }
            }
            return false;
        }

        List<PadAssignment> ClosePad(List<PadAssignment> list, int pi)
        {
            var trial = list.Select(x => x.Clone()).ToList();
            var rem = trial[pi]; trial.RemoveAt(pi);
            var order = rem.Wells.OrderBy(w => trial.Min(t => C[t.Site, w])).ThenBy(w => w).ToList();
            foreach (var w in order)
            {
                double best = double.PositiveInfinity; int bt = -1;
                for (int ti = 0; ti < trial.Count; ti++)
                {
                    var t = trial[ti];
                    if (t.Wells.Count >= K || double.IsInfinity(C[t.Site, w])) continue;
                    double d = C[t.Site, w] + PadFixed(t.Site, t.Wells.Count + 1) - PadFixed(t.Site, t.Wells.Count);
                    if (d < best || (d == best && ti < bt)) { best = d; bt = ti; }
                }
                if (bt < 0) return null;
                trial[bt].Wells.Add(w);
            }
            for (int ti = 0; ti < trial.Count; ti++) Relocate(trial, ti);
            return trial;
        }

        (List<PadAssignment>, double) LocalSearch(List<PadAssignment> list)
        {
            double cur = TotalCost(list);
            for (int it = 0; it < 200; it++)
            {
                // a) перенос скважины в другой куст
                for (int pi = 0; pi < list.Count; pi++)
                {
                    foreach (var w in list[pi].Wells.ToList())
                    {
                        var p = list[pi]; int np = p.Wells.Count;
                        if (np == 0) break;
                        double baseP = PadFixed(p.Site, np) - PadFixed(p.Site, np - 1);
                        double best = -1e-6; int bq = -1;
                        for (int qi = 0; qi < list.Count; qi++)
                        {
                            var q = list[qi];
                            if (qi == pi || q.Wells.Count == 0 || q.Wells.Count >= K || double.IsInfinity(C[q.Site, w])) continue;
                            int nq = q.Wells.Count;
                            double delta = PadFixed(q.Site, nq + 1) - PadFixed(q.Site, nq) - baseP + C[q.Site, w] - C[p.Site, w];
                            if (delta < best) { best = delta; bq = qi; }
                        }
                        if (bq >= 0) { p.Wells.Remove(w); list[bq].Wells.Add(w); }
                    }
                }
                list = list.Where(x => x.Wells.Count > 0).ToList();

                // b) обмен скважинами
                for (int pi = 0; pi < list.Count; pi++)
                    for (int qi = pi + 1; qi < list.Count; qi++)
                    {
                        var p = list[pi]; var q = list[qi];
                        foreach (var a in p.Wells.ToList())
                            foreach (var b in q.Wells.ToList())
                            {
                                if (!p.Wells.Contains(a) || !q.Wells.Contains(b)) continue;
                                double dl = C[q.Site, a] + C[p.Site, b] - C[p.Site, a] - C[q.Site, b];
                                if (!double.IsNaN(dl) && !double.IsInfinity(dl) && dl < -1e-6)
                                {
                                    p.Wells.Remove(a); q.Wells.Remove(b); p.Wells.Add(b); q.Wells.Add(a);
                                }
                            }
                    }

                // c) перенос площадок
                for (int pi = 0; pi < list.Count; pi++) Relocate(list, pi);

                // d) закрытие куста
                if (list.Count > 1)
                {
                    double curTot = TotalCost(list);
                    for (int pi = 0; pi < list.Count; pi++)
                    {
                        var trial = ClosePad(list, pi);
                        if (trial != null && TotalCost(trial) < curTot - 1e-6) { list = trial; break; }
                    }
                }

                double nc = TotalCost(list);
                if (nc > cur - 1e-6) { cur = nc; break; }
                cur = nc;
            }
            return (list, cur);
        }

        var (result, total) = LocalSearch(pads);

        // 3) принудительное сокращение числа кустов
        while (maxPads.HasValue && result.Count > maxPads.Value)
        {
            List<PadAssignment> bestTrial = null; double bestCost = double.PositiveInfinity;
            for (int pi = 0; pi < result.Count; pi++)
            {
                var trial = ClosePad(result, pi);
                if (trial == null) continue;
                double c = TotalCost(trial);
                if (c < bestCost) { bestCost = c; bestTrial = trial; }
            }
            if (bestTrial == null) return (null, double.PositiveInfinity);
            (result, total) = LocalSearch(bestTrial);
        }
        return (result, total);
    }
}
