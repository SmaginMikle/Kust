using System.Globalization;
using System.Reflection;
using System.Text;
using System.Text.Json;

namespace WellPadOptimizer;

/// <summary>Табличные результаты (CSV для Excel, консоль), JSON для визуализации и HTML-отчёт с 3D.</summary>
public static class Reporting
{
    private static readonly CultureInfo Ru = CultureInfo.GetCultureInfo("ru-RU");
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private static readonly UTF8Encoding Utf8Bom = new(true);

    private static string F(double? v, int d = 1) => v.HasValue && !double.IsInfinity(v.Value) && !double.IsNaN(v.Value) ? v.Value.ToString("F" + d, Ru) : "";
    private static string Esc(string s) => s.Contains(';') || s.Contains('"') ? "\"" + s.Replace("\"", "\"\"") + "\"" : s;

    // ------------------------------------------------------------------ консоль

    public static void PrintSummary(IEnumerable<VariantResult> variants, VariantResult best)
    {
        Console.WriteLine();
        Console.WriteLine($"{"Вариант",-52} {"Кустов",6} {"Σ MD, км",9} {"Бурение",10} {"Кусты",9} {"Отсыпка",9} {"ИТОГО",10} {"Наруш.",6}");
        Console.WriteLine($"{"",-52} {"",6} {"",9} {"млн руб",10} {"млн руб",9} {"млн руб",9} {"млн руб",10} {"",6}");
        Console.WriteLine(new string('─', 117));
        foreach (var v in variants)
        {
            string mark = ReferenceEquals(v, best) ? " ◄ рекомендуемый" : v.SameAs != null ? "  (= " + v.SameAs + ")" : "";
            string name = v.Name.Length > 52 ? v.Name[..52] : v.Name;
            Console.WriteLine($"{name,-52} {v.PadCount,6} {v.MdTotal / 1000,9:F1} {v.DrillCost / 1000,10:F1} {v.PadBuildCost / 1000,9:F1} {v.SandCost / 1000,9:F1} {v.TotalCost / 1000,10:F1} {v.Violations,6}{mark}");
        }
        Console.WriteLine();
        Console.WriteLine($"Рекомендуемый вариант: {best.Name}");
        Console.WriteLine($"{"Куст",-6} {"N",11} {"E",10} {"Скв.",5} {"Угол",5} {"Суходол",8} {"Болото",7} {"Водоём",7} {"Отсыпка м3",11} {"Строит.",9} {"Отсыпка",9}");
        foreach (var p in best.Pads)
        {
            p.TerrainFractions.TryGetValue(1, out var f1); p.TerrainFractions.TryGetValue(2, out var f2); p.TerrainFractions.TryGetValue(3, out var f3);
            Console.WriteLine($"{p.Name,-6} {p.N,11:F1} {p.E,10:F1} {p.WellCount,5} {p.Angle,5:F0} {f1,8:P0} {f2,7:P0} {f3,7:P0} {p.SandVolume,11:F0} {p.BuildCost / 1000,9:F1} {p.SandCost / 1000,9:F1}");
        }
    }

    // ------------------------------------------------------------------ CSV

    public static void WriteCsv(string dir, List<VariantResult> variants, Parameters p, List<InfeasibleVariant> infeasible)
    {
        Directory.CreateDirectory(dir);
        var sb = new StringBuilder();
        sb.AppendLine("Вариант;Сценарий;Макс. скв. на кусте;Совпадает с;Кустов;Скважин;Сумма MD, м;Бурение, тыс.руб;Строительство кустов, тыс.руб;Отсыпка, м3;Отсыпка, тыс.руб;Итого, тыс.руб;Нарушений");
        foreach (var v in variants)
            sb.AppendLine(string.Join(";", Esc(v.Name), Esc(v.Scenario ?? ""), v.K?.ToString() ?? "", Esc(v.SameAs ?? ""), v.PadCount, v.WellCount,
                F(v.MdTotal), F(v.DrillCost), F(v.PadBuildCost), F(v.SandVolume, 0), F(v.SandCost), F(v.TotalCost), v.Violations));
        foreach (var x in infeasible)
            sb.AppendLine(string.Join(";", Esc(x.Name), Esc(x.Scenario), x.K, "", "", "", "", "", "", "", "", Esc("нет решения: " + x.Reason), ""));
        File.WriteAllText(Path.Combine(dir, "summary.csv"), sb.ToString(), Utf8Bom);

        sb.Clear();
        sb.AppendLine("Вариант;Куст;N;E;Скважин;Угол ряда, °;Длина, м;Ширина, м;" +
            string.Join(";", p.Terrain.Select(t => t.Name + ", доля")) + ";Отсыпка, м3;Отсыпка, тыс.руб;Строительство, тыс.руб;Запретная зона;Скважины");
        foreach (var v in variants.Where(x => x.SameAs == null))
            foreach (var pd in v.Pads)
                sb.AppendLine(string.Join(";", Esc(v.Name), pd.Name, F(pd.N, 2), F(pd.E, 2), pd.WellCount, F(pd.Angle, 0), F(pd.Length), F(pd.Width),
                    string.Join(";", p.Terrain.Select(t => F(pd.TerrainFractions.TryGetValue(t.Code, out var f) ? f : 0, 3))),
                    F(pd.SandVolume, 0), F(pd.SandCost), F(pd.BuildCost), pd.Forbidden ? "да" : "нет", string.Join(" ", pd.Wells)));
        File.WriteAllText(Path.Combine(dir, "pads.csv"), sb.ToString(), Utf8Bom);

        sb.Clear();
        sb.AppendLine("Вариант;Скважина;UID;Куст;Слот;Устье N;Устье E;T1 N;T1 E;T1 TVD;T3 N;T3 E;T3 TVD;Отход до T1, м;Разворот, °;" +
            "MD до T1, м;MD до T3, м;ГС, м;Тангенс, м;Зенит тангенса, °;Азимут тангенса, °;Макс. ПИ, °/10м;Нагрузка на крюке, тс;Момент, кН·м;" +
            "Мин. сближение, м;ННС, тыс.руб;ГС, тыс.руб;Бурение, тыс.руб;Нарушения;Предупреждения");
        foreach (var v in variants.Where(x => x.SameAs == null))
            foreach (var w in v.Wells)
                sb.AppendLine(string.Join(";", Esc(v.Name), w.Well, w.Uid, w.Pad, w.Slot, F(w.WhN, 2), F(w.WhE, 2),
                    F(w.T1.N, 2), F(w.T1.E, 2), F(w.T1.Z, 2), F(w.T3.N, 2), F(w.T3.E, 2), F(w.T3.Z, 2), F(w.Offset), F(w.AzimuthTurn),
                    F(w.MdT1), F(w.MdT3), F(w.Horizontal), F(w.TangentLength), F(w.TangentInc, 2), F(w.TangentAz, 2), F(w.MaxDls, 2),
                    F(w.HookLoad), F(w.Torque), F(w.MinSeparation), F(w.CostNNS), F(w.CostGS), F(w.Cost),
                    string.Join(", ", w.Violations), string.Join(", ", w.Warnings)));
        File.WriteAllText(Path.Combine(dir, "wells.csv"), sb.ToString(), Utf8Bom);

        // инклинометрия по уникальным вариантам
        var trDir = Path.Combine(dir, "trajectories");
        Directory.CreateDirectory(trDir);
        int vi = 0;
        foreach (var v in variants)
        {
            vi++;
            if (v.SameAs != null) continue;
            sb.Clear();
            sb.AppendLine("Скважина;Куст;MD, м;Зенит, °;Азимут, °;N;E;TVD, м;ПИ, °/10м");
            foreach (var w in v.Wells)
                foreach (var s in w.Survey)
                    sb.AppendLine(string.Join(";", w.Well, w.Pad, F(s.MD, 2), F(s.Inc, 3), F(s.Az, 3), F(s.N, 2), F(s.E, 2), F(s.TVD, 2), F(s.Dls, 3)));
            File.WriteAllText(Path.Combine(trDir, $"variant_{vi:00}.csv"), sb.ToString(), Utf8Bom);
        }
    }

    // ------------------------------------------------------------------ JSON + HTML

    public static string BuildJson(List<VariantResult> variants, VariantResult best, InputData data, TerrainMap map, Parameters p, string inputName, string mapName, List<InfeasibleVariant> infeasible)
    {
        using var ms = new MemoryStream();
        using (var j = new Utf8JsonWriter(ms, new JsonWriterOptions { Encoder = System.Text.Encodings.Web.JavaScriptEncoder.UnsafeRelaxedJsonEscaping }))
        {
            j.WriteStartObject();
            j.WriteStartObject("meta");
            j.WriteString("input", inputName);
            j.WriteString("map", mapName);
            j.WriteString("generated", DateTime.Now.ToString("yyyy-MM-dd HH:mm", Inv));
            j.WriteString("engine", "C# WellPadOptimizer");
            j.WriteString("best", best.Name);
            j.WriteStartArray("infeasible");
            foreach (var x in infeasible)
            {
                j.WriteStartObject(); j.WriteString("name", x.Name); j.WriteString("scenario", x.Scenario);
                j.WriteNumber("K", x.K); j.WriteString("reason", x.Reason); j.WriteEndObject();
            }
            j.WriteEndArray();
            j.WriteStartObject("params");
            j.WriteNumber("MaxOffset", p.MaxOffset); j.WriteNumber("MaxMD", p.MaxMD);
            j.WriteNumber("Dls1", p.Dls1); j.WriteNumber("Dls2", p.Dls2); j.WriteNumber("VerticalLength", p.VerticalLength);
            j.WriteBoolean("UseAzimuthTurn", p.UseAzimuthTurn); j.WriteNumber("MaxAzimuthTurn", p.MaxAzimuthTurn);
            j.WriteNumber("MaxHookLoad", p.MaxHookLoad); j.WriteNumber("MaxTorque", p.MaxTorque);
            j.WriteNumber("CostNNS", p.CostNNS); j.WriteNumber("CostGS", p.CostGS); j.WriteNumber("SandCostPer1000m3", p.SandCostPer1000m3);
            j.WriteNumber("MinWellsPerPad", p.MinWellsPerPad); j.WriteNumber("StartMaxWellsPerPad", p.StartMaxWellsPerPad); j.WriteNumber("EndMaxWellsPerPad", p.EndMaxWellsPerPad);
            j.WriteEndObject();
            j.WriteEndObject();

            j.WriteNumber("forbiddenCode", p.ForbiddenCode);
            j.WriteStartArray("terrainClasses");
            foreach (var t in p.Terrain) { j.WriteStartObject(); j.WriteNumber("code", t.Code); j.WriteString("name", t.Name); j.WriteNumber("fill", t.FillHeight); j.WriteEndObject(); }
            j.WriteEndArray();

            j.WriteStartArray("regions");
            foreach (var (name, e0, n0, e1, n1) in Regions(data, p))
            {
                var grid = map.Crop(e0, n0, e1, n1, out double ce0, out double cn1, out double cell);
                j.WriteStartObject();
                j.WriteString("name", name);
                j.WriteStartArray("bounds"); j.WriteNumberValue(e0); j.WriteNumberValue(n0); j.WriteNumberValue(e1); j.WriteNumberValue(n1); j.WriteEndArray();
                j.WriteNumber("gridE0", ce0); j.WriteNumber("gridN1", cn1); j.WriteNumber("cell", cell);
                j.WriteStartArray("grid");
                foreach (var row in grid) { j.WriteStartArray(); foreach (var c in row) j.WriteNumberValue(c); j.WriteEndArray(); }
                j.WriteEndArray();
                j.WriteEndObject();
            }
            j.WriteEndArray();

            j.WriteStartArray("wells");
            foreach (var w in data.Wells)
            {
                j.WriteStartObject(); j.WriteString("name", w.Name); j.WriteString("uid", w.Uid);
                WriteVec(j, "t1", w.T1); WriteVec(j, "t3", w.T3); j.WriteEndObject();
            }
            j.WriteEndArray();

            j.WriteStartArray("variants");
            foreach (var v in variants)
            {
                j.WriteStartObject();
                j.WriteString("name", v.Name); j.WriteString("kind", v.Kind); j.WriteString("scenario", v.Scenario ?? "");
                if (v.K.HasValue) j.WriteNumber("K", v.K.Value); else j.WriteNull("K");
                if (v.SameAs != null) j.WriteString("sameAs", v.SameAs); else j.WriteNull("sameAs");
                j.WriteNumber("nPads", v.PadCount); j.WriteNumber("nWells", v.WellCount);
                Num(j, "mdTotal", v.MdTotal); Num(j, "drillCost", v.DrillCost); Num(j, "padBuildCost", v.PadBuildCost);
                Num(j, "sandCost", v.SandCost); Num(j, "sandVolume", v.SandVolume); Num(j, "totalCost", v.TotalCost);
                j.WriteNumber("violations", v.Violations);
                if (v.OptimizerEstimate.HasValue) Num(j, "estimate", v.OptimizerEstimate.Value); else j.WriteNull("estimate");
                if (v.SameAs == null)
                {
                    j.WriteStartArray("pads");
                    foreach (var pd in v.Pads)
                    {
                        j.WriteStartObject();
                        j.WriteString("name", pd.Name); Num(j, "n", pd.N, 2); Num(j, "e", pd.E, 2); j.WriteNumber("count", pd.WellCount);
                        Num(j, "angle", pd.Angle); Num(j, "length", pd.Length); Num(j, "width", pd.Width);
                        j.WriteStartObject("terrain");
                        foreach (var kv in pd.TerrainFractions) Num(j, kv.Key.ToString(Inv), kv.Value, 3);
                        j.WriteEndObject();
                        Num(j, "sandVolume", pd.SandVolume, 0); Num(j, "sandCost", pd.SandCost); Num(j, "buildCost", pd.BuildCost);
                        j.WriteBoolean("forbidden", pd.Forbidden);
                        j.WriteStartArray("wells"); foreach (var s in pd.Wells) j.WriteStringValue(s); j.WriteEndArray();
                        j.WriteEndObject();
                    }
                    j.WriteEndArray();
                    j.WriteStartArray("wells");
                    foreach (var w in v.Wells)
                    {
                        j.WriteStartObject();
                        j.WriteString("well", w.Well); j.WriteString("pad", w.Pad); j.WriteNumber("slot", w.Slot);
                        Num(j, "whN", w.WhN, 2); Num(j, "whE", w.WhE, 2); Num(j, "offset", w.Offset); Num(j, "azTurn", w.AzimuthTurn);
                        Num(j, "hz", w.Horizontal); NumN(j, "mdT1", w.MdT1); NumN(j, "mdT3", w.MdT3); NumN(j, "tanLen", w.TangentLength);
                        NumN(j, "tanInc", w.TangentInc, 2); NumN(j, "tanAz", w.TangentAz, 2); NumN(j, "cost", w.Cost); NumN(j, "costNNS", w.CostNNS);
                        NumN(j, "costGS", w.CostGS); NumN(j, "hook", w.HookLoad); NumN(j, "torque", w.Torque); NumN(j, "maxDls", w.MaxDls, 3);
                        NumN(j, "minSep", w.MinSeparation);
                        j.WriteStartArray("viol"); foreach (var s in w.Violations) j.WriteStringValue(s); j.WriteEndArray();
                        j.WriteStartArray("warn"); foreach (var s in w.Warnings) j.WriteStringValue(s); j.WriteEndArray();
                        j.WriteStartArray("tr");
                        for (int i = 0; i < w.Survey.Count; i++)
                        {
                            var s = w.Survey[i];
                            bool key = i % 5 == 0 || i == w.Survey.Count - 1 || (w.MdT1.HasValue && Math.Abs(s.MD - w.MdT1.Value) < 1e-6);
                            if (!key) continue;
                            j.WriteStartArray();
                            j.WriteNumberValue(Math.Round(s.N, 1)); j.WriteNumberValue(Math.Round(s.E, 1)); j.WriteNumberValue(Math.Round(s.TVD, 1));
                            j.WriteNumberValue(Math.Round(s.MD, 1));
                            j.WriteEndArray();
                        }
                        j.WriteEndArray();
                        j.WriteEndObject();
                    }
                    j.WriteEndArray();
                }
                j.WriteEndObject();
            }
            j.WriteEndArray();
            j.WriteEndObject();
        }
        return Encoding.UTF8.GetString(ms.ToArray());
    }

    private static void WriteVec(Utf8JsonWriter j, string name, Vec3 v)
    {
        j.WriteStartArray(name); j.WriteNumberValue(Math.Round(v.N, 2)); j.WriteNumberValue(Math.Round(v.E, 2)); j.WriteNumberValue(Math.Round(v.Z, 2)); j.WriteEndArray();
    }
    private static void Num(Utf8JsonWriter j, string name, double v, int d = 1) => j.WriteNumber(name, Math.Round(v, d));
    private static void NumN(Utf8JsonWriter j, string name, double? v, int d = 1)
    {
        if (v.HasValue && !double.IsInfinity(v.Value) && !double.IsNaN(v.Value)) j.WriteNumber(name, Math.Round(v.Value, d)); else j.WriteNull(name);
    }

    /// <summary>Участки: кластеры целей, удалённые друг от друга более чем на 5 км (односвязная кластеризация).</summary>
    public static List<(string name, double e0, double n0, double e1, double n1)> Regions(InputData data, Parameters p)
    {
        var w = data.Wells; int n = w.Count;
        var lab = Enumerable.Range(0, n).ToArray();
        int Find(int x) { while (lab[x] != x) x = lab[x] = lab[lab[x]]; return x; }
        for (int i = 0; i < n; i++)
            for (int k = i + 1; k < n; k++)
                if (Geo.Hypot(w[i].T1.N - w[k].T1.N, w[i].T1.E - w[k].T1.E) < 5000) lab[Find(i)] = Find(k);
        var groups = Enumerable.Range(0, n).GroupBy(Find).Select(g => g.ToList())
            .OrderBy(g => g.Average(i => w[i].T1.N)).ToList();
        var res = new List<(string, double, double, double, double)>();
        double m = p.MaxOffset + 300;
        int no = 1;
        foreach (var g in groups)
        {
            double e0 = g.Min(i => Math.Min(w[i].T1.E, w[i].T3.E)) - m, e1 = g.Max(i => Math.Max(w[i].T1.E, w[i].T3.E)) + m;
            double n0 = g.Min(i => Math.Min(w[i].T1.N, w[i].T3.N)) - m, n1 = g.Max(i => Math.Max(w[i].T1.N, w[i].T3.N)) + m;
            res.Add(($"Участок {no++} ({g.Count} скв.)", Math.Floor(e0 / 100) * 100, Math.Floor(n0 / 100) * 100, Math.Ceiling(e1 / 100) * 100, Math.Ceiling(n1 / 100) * 100));
        }
        return res;
    }

    public static void WriteHtml(string path, string json)
    {
        string tpl;
        var asm = Assembly.GetExecutingAssembly();
        using (var s = asm.GetManifestResourceStream("report_template.html"))
        {
            if (s == null)
            {
                var local = Path.Combine(AppContext.BaseDirectory, "report_template.html");
                tpl = File.ReadAllText(local);
            }
            else using (var r = new StreamReader(s)) tpl = r.ReadToEnd();
        }
        // защита от закрывающего тега внутри JSON
        json = json.Replace("</", "<\\/");
        string body = tpl.Replace("/*__REPORT_DATA__*/null", json);
        // шаблон — фрагмент страницы; для автономного файла добавляем полноценный каркас документа
        string html = "<!doctype html>\n<html lang=\"ru\">\n<head>\n<meta charset=\"utf-8\">\n" +
                      "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n" +
                      "<style>body{margin:0}[hidden]{display:none!important}</style>\n</head>\n<body>\n" + body + "\n</body>\n</html>\n";
        File.WriteAllText(path, html, new UTF8Encoding(false));
    }
}
