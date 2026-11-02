namespace WellPadOptimizer;

/// <summary>
/// Проектный профиль: вертикаль (VerticalLength) → первый набор (R1) → тангенциальный участок →
/// второй набор/разворот (R2) с выходом в T1 по направлению T1→T3 → горизонтальный участок T1→T3.
/// Для целей без горизонтального участка (T1 = T3) — J-профиль: вертикаль → набор → тангенс до T1.
/// </summary>
public sealed class Profile
{
    public Vec3 Wellhead;          // устье (Z = 0)
    public Vec3 Kop;               // точка начала набора
    public Vec3 TangentDir;        // направление тангенциального участка (ед. вектор)
    public Vec3? TargetDir;        // направление T1→T3 (null для ННС без ГС)
    public double Theta1, Theta2;  // углы искривления на участках набора, рад
    public double TangentLength;
    public double MdT1, MdT3, Horizontal;
    public double TangentInclination => Math.Acos(Geo.Clamp1(TangentDir.Z)) / Geo.Deg;
    public double TangentAzimuth => Geo.Azimuth(TangentDir.N, TangentDir.E);
}

/// <summary>Точка инклинометрии.</summary>
public readonly record struct SurveyPoint(double MD, double Inc, double Az, double N, double E, double TVD, double Dls)
{
    public Vec3 Pos => new(N, E, TVD);
}

public enum CheckResult { Ok, Offset, AzimuthTurn, Profile, MD }

public static class TrajectoryCalculator
{
    /// <summary>
    /// Решение профиля с заданным выходом в T1. Неизвестное — направление тангенса d,
    /// находится итерацией неподвижной точки: d = норм(P2(d) − P1(d)), где
    /// P1 = KOP + R1·tg(θ1/2)·(v + d) — конец первого набора,
    /// P2 = T1 − R2·tg(θ2/2)·(d + t) — начало второго набора (метод минимальной кривизны).
    /// </summary>
    public static Profile Solve(Parameters p, Vec3 wellhead, Vec3 t1, Vec3 t3)
    {
        double r1 = p.R1, r2 = p.R2;
        var v = Vec3.Down;
        var kop = new Vec3(wellhead.N, wellhead.E, wellhead.Z + p.VerticalLength);
        double hz = (t3 - t1).Length;
        Vec3? t = hz > 1e-6 ? (t3 - t1).Unit() : null;
        var d = (t1 - kop).Unit();
        bool ok = false;
        for (int it = 0; it < 300; it++)
        {
            Ends(d, out var p1, out var p2, out _, out _, out bool bad);
            if (bad) return null;
            var w = p2 - p1; double lw = w.Length;
            if (lw < 1e-6) return null;
            var dn = w * (1.0 / lw);
            if ((dn - d).Length < 1e-10) { d = dn; ok = true; break; }
            d = (0.5 * d + 0.5 * dn).Unit();
        }
        if (!ok) return null;
        Ends(d, out var e1, out var e2, out double th1, out double th2, out bool bad2);
        if (bad2) return null;
        double L = (e2 - e1).Dot(d);
        if (L < -1e-6 || d.Z < -1e-9) return null;
        double mdT1 = p.VerticalLength + r1 * th1 + Math.Max(L, 0) + r2 * th2;
        return new Profile
        {
            Wellhead = wellhead, Kop = kop, TangentDir = d, TargetDir = t, Theta1 = th1, Theta2 = th2,
            TangentLength = Math.Max(L, 0), MdT1 = mdT1, MdT3 = mdT1 + hz, Horizontal = hz,
        };

        void Ends(Vec3 dd, out Vec3 a1, out Vec3 a2, out double q1, out double q2, out bool fail)
        {
            fail = false;
            q1 = Math.Acos(Geo.Clamp1(dd.Dot(v)));
            a1 = kop + (r1 * Math.Tan(q1 / 2)) * (v + dd);
            if (t.HasValue)
            {
                q2 = Math.Acos(Geo.Clamp1(dd.Dot(t.Value)));
                if (q2 > Math.PI - 1e-3) { fail = true; a2 = t1; return; }
                a2 = t1 - (r2 * Math.Tan(q2 / 2)) * (dd + t.Value);
            }
            else { q2 = 0; a2 = t1; }
        }
    }

    /// <summary>Угол разворота: между азимутом «устье→T1» и азимутом горизонтального участка T1→T3.</summary>
    public static double AzimuthTurn(Vec3 wh, Vec3 t1, Vec3 t3)
    {
        if ((t3 - t1).HorizontalLength < 1e-6) return 0;
        double a1 = Geo.Azimuth(t1.N - wh.N, t1.E - wh.E);
        double a2 = Geo.Azimuth(t3.N - t1.N, t3.E - t1.E);
        return Geo.AngleDiff(a1, a2);
    }

    /// <summary>Полная проверка технологических ограничений для пары устье–цель.</summary>
    public static (Profile prof, CheckResult res) Check(Parameters p, Vec3 wh, WellTarget w)
    {
        var lim = p.LimitsToBottomHole ? w.T3 : w.T1;
        if (Geo.Hypot(lim.N - wh.N, lim.E - wh.E) > p.MaxOffset + 1e-6) return (null, CheckResult.Offset);
        if (p.UseAzimuthTurn && AzimuthTurn(wh, w.T1, w.T3) > p.MaxAzimuthTurn + 1e-6) return (null, CheckResult.AzimuthTurn);
        var prof = Solve(p, wh, w.T1, w.T3);
        if (prof == null) return (null, CheckResult.Profile);
        double md = p.LimitsToBottomHole ? prof.MdT3 : prof.MdT1;
        if (md > p.MaxMD + 1e-6) return (null, CheckResult.MD);
        return (prof, CheckResult.Ok);
    }

    /// <summary>Стоимость бурения скважины, тыс. руб.: ННС до T1 + ГС T1–T3.</summary>
    public static double DrillCost(Parameters p, Profile prof) => p.CostNNS * prof.MdT1 + p.CostGS * prof.Horizontal;

    // ------------------------------------------------------------------ инклинометрия

    public static List<SurveyPoint> Survey(Parameters p, Profile pr, double step = 0)
    {
        if (step <= 0) step = p.TrajectoryStep;
        var raw = new List<(double md, Vec3 pos, Vec3 tan)>();
        var v = Vec3.Down;
        var cur = pr.Wellhead; double md = 0;
        raw.Add((0, cur, v));
        Line(raw, ref cur, ref md, v, p.VerticalLength, step);
        Arc(raw, ref cur, ref md, v, pr.TangentDir, p.R1, step);
        Line(raw, ref cur, ref md, pr.TangentDir, pr.TangentLength, step);
        if (pr.TargetDir.HasValue)
        {
            Arc(raw, ref cur, ref md, pr.TangentDir, pr.TargetDir.Value, p.R2, step);
            Line(raw, ref cur, ref md, pr.TargetDir.Value, pr.Horizontal, step);
        }
        var res = new List<SurveyPoint>(raw.Count);
        for (int i = 0; i < raw.Count; i++)
        {
            var (m, pos, tg) = raw[i];
            double inc = Math.Acos(Geo.Clamp1(tg.Z)) / Geo.Deg;
            double az = inc > 1e-6 ? Geo.Azimuth(tg.N, tg.E) : 0;
            double dls = 0;
            if (i > 0 && m - raw[i - 1].md > 1e-6)
                dls = Math.Acos(Geo.Clamp1(raw[i - 1].tan.Dot(tg))) / Geo.Deg / (m - raw[i - 1].md) * p.DlsBaseLength;
            res.Add(new SurveyPoint(m, inc, az, pos.N, pos.E, pos.Z, dls));
        }
        return res;
    }

    private static void Line(List<(double, Vec3, Vec3)> o, ref Vec3 cur, ref double md, Vec3 d, double len, double step)
    {
        if (len <= 1e-9) return;
        int k = 1;
        for (; k * step < len - 1e-6; k++) o.Add((md + k * step, cur + (k * step) * d, d));
        cur = cur + len * d; md += len;
        o.Add((md, cur, d));
    }

    private static void Arc(List<(double, Vec3, Vec3)> o, ref Vec3 cur, ref double md, Vec3 a, Vec3 b, double r, double step)
    {
        double th = Math.Acos(Geo.Clamp1(a.Dot(b)));
        if (th < 1e-9) return;
        var n = (b - a.Dot(b) * a).Unit();
        double larc = r * th;
        var p0 = cur;
        for (int k = 1; k * step < larc - 1e-6; k++)
        {
            double ph = k * step / r;
            o.Add((md + k * step, p0 + r * (Math.Sin(ph) * a + (1 - Math.Cos(ph)) * n), Math.Cos(ph) * a + Math.Sin(ph) * n));
        }
        cur = p0 + r * (Math.Sin(th) * a + (1 - Math.Cos(th)) * n);
        md += larc;
        o.Add((md, cur, b));
    }

    // ------------------------------------------------------------------ нагрузки

    /// <summary>
    /// «Мягкая» модель колонны (Johancsik): нагрузка на крюке при подъёме, тс, и момент
    /// при вращении над забоем, кН·м. Интегрирование от забоя (T3) к устью.
    /// </summary>
    public static (double hookLoadTf, double torqueKNm) TorqueDrag(Parameters p, List<SurveyPoint> sv, double mdT1)
    {
        double bf = 1 - p.MudDensity / p.SteelDensity;
        double w = p.PipeWeightAir * bf;
        double tPull = 0, tRot = 0, tq = 0;
        for (int i = sv.Count - 1; i > 0; i--)
        {
            var b = sv[i]; var a = sv[i - 1];
            double dL = b.MD - a.MD; if (dL <= 0) continue;
            double i1 = a.Inc * Geo.Deg, i2 = b.Inc * Geo.Deg;
            double ia = (i1 + i2) / 2, di = i2 - i1;
            double da = (b.Az - a.Az) * Geo.Deg;
            da = ((da + Math.PI) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI) - Math.PI;
            if (a.Inc < 1e-6 || b.Inc < 1e-6) da = 0;
            double mu = a.MD >= mdT1 ? p.FrictionOpenHole : p.FrictionCased;
            double nPull = Geo.Hypot(tPull * da * Math.Sin(ia), tPull * di + w * dL * Math.Sin(ia));
            tPull += w * dL * Math.Cos(ia) + mu * nPull;
            double nRot = Geo.Hypot(tRot * da * Math.Sin(ia), tRot * di + w * dL * Math.Sin(ia));
            tRot += w * dL * Math.Cos(ia);
            tq += mu * nRot * p.PipeRadius;
        }
        return (tPull / 9.80665, tq);
    }
}
