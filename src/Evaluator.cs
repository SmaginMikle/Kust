namespace WellPadOptimizer;

public sealed class PadResult
{
    public int Id { get; set; }
    public string Name { get; set; } = "";
    public double N { get; set; }
    public double E { get; set; }
    public int WellCount { get; set; }
    public double Angle { get; set; }
    public double Length { get; set; }
    public double Width { get; set; }
    public Dictionary<int, double> TerrainFractions { get; set; } = new();
    public double SandVolume { get; set; }
    public double SandCost { get; set; }
    public double BuildCost { get; set; }
    public bool Forbidden { get; set; }
    public List<string> Wells { get; set; } = new();
}

public sealed class WellResult
{
    public string Well { get; set; } = "";
    public string Uid { get; set; } = "";
    public string Pad { get; set; } = "";
    public int Slot { get; set; }
    public double WhN { get; set; }
    public double WhE { get; set; }
    public double Offset { get; set; }
    public double AzimuthTurn { get; set; }
    public Vec3 T1 { get; set; }
    public Vec3 T3 { get; set; }
    public double Horizontal { get; set; }
    public double? MdT1 { get; set; }
    public double? MdT3 { get; set; }
    public double? TangentLength { get; set; }
    public double? TangentInc { get; set; }
    public double? TangentAz { get; set; }
    public double? Cost { get; set; }
    public double? CostNNS { get; set; }
    public double? CostGS { get; set; }
    public double? HookLoad { get; set; }
    public double? Torque { get; set; }
    public double? MaxDls { get; set; }
    public double? MinSeparation { get; set; }
    public List<string> Violations { get; set; } = new();
    public List<string> Warnings { get; set; } = new();
    public List<SurveyPoint> Survey { get; set; } = new();
}

public sealed class VariantResult
{
    public string Name { get; set; } = "";
    public string Kind { get; set; } = "";        // input | opt
    public string Scenario { get; set; } = "";
    public int? K { get; set; }
    public string SameAs { get; set; }
    public double? OptimizerEstimate { get; set; }
    public List<PadResult> Pads { get; set; } = new();
    public List<WellResult> Wells { get; set; } = new();
    public int PadCount => Pads.Count;
    public int WellCount => Wells.Count;
    public double MdTotal => Wells.Sum(w => w.MdT3 ?? 0);
    public double DrillCost => Wells.Sum(w => w.Cost ?? 0);
    public double PadBuildCost => Pads.Sum(p => p.BuildCost);
    public double SandCost => Pads.Sum(p => p.SandCost);
    public double SandVolume => Pads.Sum(p => p.SandVolume);
    public double TotalCost => DrillCost + PadBuildCost + SandCost;
    public int Violations => Wells.Count(w => w.Violations.Count > 0) + Pads.Count(p => p.Forbidden);

    public VariantResult AliasOf(string name, int? k) => new()
    {
        Name = name, Kind = Kind, Scenario = Scenario, K = k, SameAs = Name,
        OptimizerEstimate = OptimizerEstimate, Pads = Pads, Wells = Wells,
    };
}

/// <summary>Расчёт варианта: устья по слотам, траектории, инклинометрия, нагрузки, антиколлизия, затраты.</summary>
public sealed class Evaluator
{
    private readonly Parameters _p;
    private readonly InputData _data;
    private readonly TerrainMap _map;
    private readonly List<(double s, double w)> _samples;

    public Evaluator(Parameters p, InputData data, TerrainMap map)
    {
        _p = p; _data = data; _map = map;
        _samples = Optimizer.FootprintSamples(p, p.MaxPadCapacity);
    }

    /// <param name="pads">центр куста (N, E) и индексы скважин</param>
    /// <param name="fixedWellheads">заданные устья (для исходных вариантов) или null</param>
    public VariantResult Evaluate(string name, string kind, string scenario, int? k,
        List<((double n, double e) center, List<int> wells)> pads, Dictionary<int, (double n, double e)> fixedWellheads = null)
    {
        var res = new VariantResult { Name = name, Kind = kind, Scenario = scenario, K = k };
        int padNo = _p.StartPadNumber;
        foreach (var (center, wells) in pads)
        {
            int n = wells.Count;
            double angle; (bool forbidden, Dictionary<int, double> frac, double volume, double cost) fp;
            var slots = new List<(int well, double n, double e)>();
            if (fixedWellheads != null)
            {
                angle = PrincipalAngle(wells.Select(w => fixedWellheads[w]).ToList());
                fp = Optimizer.EvaluateFootprint(_p, _map, _samples, center.n, center.e, n, angle);
                foreach (var w in wells) slots.Add((w, fixedWellheads[w].n, fixedWellheads[w].e));
            }
            else
            {
                (angle, fp) = ChooseOrientation(center, wells);
                double uN = Math.Cos(angle * Geo.Deg), uE = Math.Sin(angle * Geo.Deg);
                var offs = _p.SlotOffsets(n);
                var order = wells.OrderBy(w => (_data.Wells[w].T1.N - center.n) * uN + (_data.Wells[w].T1.E - center.e) * uE)
                                 .ThenBy(w => w).ToList();
                order = ImproveSlots(order, center, uN, uE, offs);
                for (int i = 0; i < n; i++) slots.Add((order[i], center.n + offs[i] * uN, center.e + offs[i] * uE));
            }

            var pad = new PadResult
            {
                Id = padNo, Name = _p.PadPrefix + padNo, N = center.n, E = center.e, WellCount = n, Angle = angle,
                Length = _p.PadLength(n), Width = _p.PadWidth, TerrainFractions = fp.frac,
                SandVolume = fp.volume, SandCost = fp.cost, BuildCost = _p.PadBuildCost(n), Forbidden = fp.forbidden,
            };
            padNo++;
            int slotNo = 1;
            foreach (var (w, sn, se) in slots)
            {
                var tgt = _data.Wells[w];
                var wh = new Vec3(sn, se, 0);
                var wr = new WellResult
                {
                    Well = tgt.Name, Uid = tgt.Uid, Pad = pad.Name, Slot = slotNo++, WhN = sn, WhE = se,
                    Offset = Geo.Hypot(tgt.T1.N - sn, tgt.T1.E - se),
                    AzimuthTurn = TrajectoryCalculator.AzimuthTurn(wh, tgt.T1, tgt.T3),
                    T1 = tgt.T1, T3 = tgt.T3, Horizontal = tgt.HorizontalLength,
                };
                var lim = _p.LimitsToBottomHole ? tgt.T3 : tgt.T1;
                if (Geo.Hypot(lim.N - sn, lim.E - se) > _p.MaxOffset + 1e-6) wr.Violations.Add("отход");
                if (_p.UseAzimuthTurn && wr.AzimuthTurn > _p.MaxAzimuthTurn + 1e-6) wr.Violations.Add("разворот");
                var prof = TrajectoryCalculator.Solve(_p, wh, tgt.T1, tgt.T3);
                if (prof == null) wr.Violations.Add("профиль");
                else
                {
                    var sv = TrajectoryCalculator.Survey(_p, prof);
                    var (hook, tq) = TrajectoryCalculator.TorqueDrag(_p, sv, prof.MdT1);
                    double md = _p.LimitsToBottomHole ? prof.MdT3 : prof.MdT1;
                    if (md > _p.MaxMD + 1e-6) wr.Violations.Add("MD");
                    if (hook > _p.MaxHookLoad) wr.Violations.Add("крюк");
                    if (tq > _p.MaxTorque) wr.Violations.Add("момент");
                    wr.MdT1 = prof.MdT1; wr.MdT3 = prof.MdT3; wr.TangentLength = prof.TangentLength;
                    wr.TangentInc = prof.TangentInclination; wr.TangentAz = prof.TangentAzimuth;
                    wr.CostNNS = _p.CostNNS * prof.MdT1; wr.CostGS = _p.CostGS * prof.Horizontal;
                    wr.Cost = wr.CostNNS + wr.CostGS;
                    wr.HookLoad = hook; wr.Torque = tq; wr.MaxDls = sv.Max(x => x.Dls);
                    wr.Survey = sv;
                }
                pad.Wells.Add(wr.Well);
                res.Wells.Add(wr);
            }
            res.Pads.Add(pad);
        }
        AntiCollision(res);
        return res;
    }

    /// <summary>Ориентация площадки: без запретных зон, минимум отсыпки, затем — поперёк среднего азимута ГС.</summary>
    private (double angle, (bool, Dictionary<int, double>, double, double) fp) ChooseOrientation((double n, double e) c, List<int> wells)
    {
        double? meanAz = null;
        double sx = 0, sy = 0; int cnt = 0;
        foreach (var w in wells)
        {
            var t = _data.Wells[w];
            if (!t.HasHorizontal) continue;
            double a = Geo.Azimuth(t.T3.N - t.T1.N, t.T3.E - t.T1.E) * Geo.Deg;
            sx += Math.Cos(2 * a); sy += Math.Sin(2 * a); cnt++;
        }
        if (cnt > 0) meanAz = ((Math.Atan2(sy, sx) / Geo.Deg / 2) % 180 + 180) % 180;
        (bool, double, double) bestKey = (true, double.MaxValue, double.MaxValue);
        double bestAng = 0; (bool, Dictionary<int, double>, double, double) bestFp = default;
        bool first = true;
        foreach (var a in _p.PadOrientations)
        {
            var fp = Optimizer.EvaluateFootprint(_p, _map, _samples, c.n, c.e, wells.Count, a);
            double perp = meanAz.HasValue ? Math.Abs((((a - (meanAz.Value + 90)) + 90) % 180 + 180) % 180 - 90) : 0;
            var key = (fp.forbidden, Math.Round(fp.cost), perp);
            if (first || Compare(key, bestKey) < 0) { bestKey = key; bestAng = a; bestFp = fp; first = false; }
        }
        return (bestAng, bestFp);

        static int Compare((bool f, double c, double p) x, (bool f, double c, double p) y)
        {
            if (x.f != y.f) return x.f ? 1 : -1;
            if (x.c != y.c) return x.c.CompareTo(y.c);
            return x.p.CompareTo(y.p);
        }
    }

    private static double PrincipalAngle(List<(double n, double e)> pts)
    {
        if (pts.Count < 2) return 0;
        double mn = pts.Average(p => p.n), me = pts.Average(p => p.e);
        double snn = 0, see = 0, sne = 0;
        foreach (var (n, e) in pts) { snn += (n - mn) * (n - mn); see += (e - me) * (e - me); sne += (n - mn) * (e - me); }
        double th = 0.5 * Math.Atan2(2 * sne, snn - see);          // угол главной оси от оси N
        double az = Geo.Azimuth(Math.Cos(th), Math.Sin(th));
        return az % 180;
    }

    /// <summary>
    /// Раскладка скважин по слотам: попарные перестановки, максимизирующие минимальное
    /// сближение стволов (до AntiCollisionMaxMD, по равным MD), затем минимальную стоимость бурения.
    /// </summary>
    private List<int> ImproveSlots(List<int> order, (double n, double e) c, double uN, double uE, double[] offs)
    {
        int n = order.Count;
        if (n < 2) return order;
        var cache = new Dictionary<(int, int), (Vec3[] pts, double cost)>();
        (Vec3[] pts, double cost) GetPath(int w, int k)
        {
            if (cache.TryGetValue((w, k), out var v)) return v;
            var t = _data.Wells[w];
            var wh = new Vec3(c.n + offs[k] * uN, c.e + offs[k] * uE, 0);
            var pr = TrajectoryCalculator.Solve(_p, wh, t.T1, t.T3);
            if (pr == null) v = (null, double.PositiveInfinity);
            else
            {
                var sv = TrajectoryCalculator.Survey(_p, pr);
                v = (sv.Where(x => x.MD > _p.VerticalLength && x.MD <= _p.AntiCollisionMaxMD).Select(x => x.Pos).ToArray(),
                     TrajectoryCalculator.DrillCost(_p, pr));
            }
            cache[(w, k)] = v;
            return v;
        }
        (double sep, double negCost) Score(List<int> ord)
        {
            double ms = double.PositiveInfinity, cost = 0;
            for (int i = 0; i < n; i++)
            {
                var a = GetPath(ord[i], i); cost += a.cost;
                for (int j = i + 1; j < n; j++) ms = Math.Min(ms, SameMdSeparation(a.pts, GetPath(ord[j], j).pts));
            }
            return (Math.Min(ms, _p.MinSeparation), -cost);
        }
        static bool Better((double s, double c) x, (double s, double c) y) => x.s > y.s || (x.s == y.s && x.c > y.c);

        var cur = Score(order);
        order = new List<int>(order);
        for (int pass = 0; pass < 4; pass++)
        {
            bool imp = false;
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                {
                    (order[i], order[j]) = (order[j], order[i]);
                    var sc = Score(order);
                    if (Better(sc, cur)) { cur = sc; imp = true; }
                    else (order[i], order[j]) = (order[j], order[i]);
                }
            if (!imp) break;
        }
        return order;
    }

    private static double SameMdSeparation(Vec3[] a, Vec3[] b)
    {
        if (a == null || b == null || a.Length == 0 || b.Length == 0) return double.PositiveInfinity;
        int m = Math.Min(a.Length, b.Length); double best = double.PositiveInfinity;
        for (int i = 0; i < m; i++) best = Math.Min(best, (a[i] - b[i]).Length);
        return best;
    }

    /// <summary>Минимальное расстояние между стволами одного куста (ниже вертикального участка).</summary>
    private void AntiCollision(VariantResult res)
    {
        foreach (var grp in res.Wells.GroupBy(w => w.Pad))
        {
            var list = grp.ToList();
            for (int i = 0; i < list.Count; i++)
                for (int j = i + 1; j < list.Count; j++)
                {
                    var a = list[i].Survey.Where(x => x.MD > _p.VerticalLength).Select(x => x.Pos).ToArray();
                    var b = list[j].Survey.Where(x => x.MD > _p.VerticalLength).Select(x => x.Pos).ToArray();
                    if (a.Length == 0 || b.Length == 0) continue;
                    double dm = double.PositiveInfinity;
                    foreach (var pa in a)
                        foreach (var pb in b)
                        {
                            double dn = pa.N - pb.N; if (Math.Abs(dn) >= dm) continue;
                            double de = pa.E - pb.E; if (Math.Abs(de) >= dm) continue;
                            double d = Math.Sqrt(dn * dn + de * de + (pa.Z - pb.Z) * (pa.Z - pb.Z));
                            if (d < dm) dm = d;
                        }
                    foreach (var r in new[] { list[i], list[j] })
                        r.MinSeparation = r.MinSeparation.HasValue ? Math.Min(r.MinSeparation.Value, dm) : dm;
                }
        }
        foreach (var w in res.Wells)
            if (w.MinSeparation.HasValue && w.MinSeparation.Value < _p.SeparationWarning)
                w.Warnings.Add($"сближение {w.MinSeparation.Value:F1} м");
    }
}
