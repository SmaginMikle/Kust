using System.Diagnostics;
using System.Text;

namespace WellPadOptimizer;

/// <summary>
/// Кустование скважин и расчёт траекторий.
/// Запуск:  WellPadOptimizer --input data/input_mass.json [--map data/terrain_map.asc] [--settings settings.json] [--out results]
/// </summary>
/// <summary>Вариант, для которого нет допустимого разбиения скважин на кусты.</summary>
public sealed record InfeasibleVariant(string Name, string Scenario, int K, string Reason);

public static class Program
{
    public static int Main(string[] args)
    {
        Console.OutputEncoding = Encoding.UTF8;
        string input = Arg(args, "--input") ?? Path.Combine(AppContext.BaseDirectory, "data", "input_mass.json");
        string mapPath = Arg(args, "--map") ?? Path.Combine(AppContext.BaseDirectory, "data", "terrain_map.asc");
        string settings = Arg(args, "--settings") ?? Path.Combine(AppContext.BaseDirectory, "settings.json");
        string outDir = Arg(args, "--out") ?? "results";
        Directory.CreateDirectory(outDir);
        var sw = Stopwatch.StartNew();

        try
        {
            var p = Parameters.Load(settings);
            var data = InputData.Load(input);
            Console.WriteLine($"Целей: {data.Wells.Count}, исходных вариантов: {data.Variants.Count}, фактических скважин: {data.RealWells.Count}");
            Console.WriteLine($"R1 = {p.R1:F0} м, R2 = {p.R2:F0} м (темпы {p.Dls1} и {p.Dls2} °/{p.DlsBaseLength} м)");

            TerrainMap map;
            if (File.Exists(mapPath)) { map = TerrainMap.LoadAsc(mapPath); }
            else
            {
                Console.WriteLine("Карта желательности не задана — построена синтетическая карта (сохранена в terrain_map_synthetic.asc)");
                map = TerrainMap.Synthetic(data.Wells);
                mapPath = Path.Combine(outDir, "terrain_map_synthetic.asc");
                map.SaveAsc(mapPath);
            }
            map.ForbiddenCode = p.ForbiddenCode;

            var opt = new Optimizer(p, data, map);
            opt.BuildSites();
            opt.BuildSandTable();
            Console.WriteLine($"Кандидатных площадок: {opt.NSites}; таблица отсыпки готова ({sw.Elapsed.TotalSeconds:F1} с)");

            var results = new List<VariantResult>();
            var evalMain = new Evaluator(p, data, map);

            // 1) исходные варианты из файла — расчёт по той же модели затрат
            foreach (var iv in data.Variants)
            {
                var pads = new List<((double, double), List<int>)>();
                var wh = new Dictionary<int, (double n, double e)>();
                foreach (var kv in iv.Pads)
                {
                    foreach (var x in kv.Value) wh[x.well] = (x.n, x.e);
                    pads.Add(((kv.Value.Average(x => x.n), kv.Value.Average(x => x.e)), kv.Value.Select(x => x.well).ToList()));
                }
                var r = evalMain.Evaluate($"Исходный вариант {iv.Id}", "input", "Исходные данные", null, pads, wh);
                results.Add(r);
                Console.WriteLine($"{r.Name}: {r.TotalCost / 1000:F1} млн руб, нарушений {r.Violations}");
            }

            // 2) перебор максимального числа скважин на кусте
            var scenarios = new List<(string name, Parameters prm)>
            {
                (p.UseAzimuthTurn ? $"Разворот ≤{p.MaxAzimuthTurn:0}°" : "Без ограничения разворота", p)
            };
            if (p.UseAzimuthTurn && p.RunScenarioWithoutAzimuthTurn)
            {
                var p2 = p.Clone(); p2.UseAzimuthTurn = false;
                scenarios.Add(("Без ограничения разворота", p2));
            }
            var infeasible = new List<InfeasibleVariant>();
            var sigs = new List<(List<(int, int[])> sig, VariantResult res)>();
            foreach (var (scName, prm) in scenarios)
            {
                var o = new Optimizer(prm, data, map) { SiteN = opt.SiteN, SiteE = opt.SiteE, Sand = opt.Sand };
                o.BuildCostMatrix();
                Console.WriteLine($"[{scName}] пар площадка–скважина: {o.FeasiblePairs}, из них без нарушений: {o.CleanPairs}");
                var ev = new Evaluator(prm, data, map);

                VariantResult AddVariant(string name, string scenario, int K, List<PadAssignment> pads, double est)
                {
                    var sig = pads.Select(x => (x.Site, x.Wells.OrderBy(w => w).ToArray())).OrderBy(x => x.Site).ToList();
                    foreach (var (s0, r0) in sigs)
                        if (SameSig(sig, s0)) { var al = r0.AliasOf(name, K); al.Scenario = scenario; results.Add(al); return al; }
                    var r = ev.Evaluate(name, "opt", scenario, K, ToPads(o, pads));
                    r.OptimizerEstimate = est;
                    results.Add(r); sigs.Add((sig, r));
                    Console.WriteLine($"{name}: кустов {r.PadCount} [{string.Join(",", r.Pads.Select(x => x.WellCount))}], " +
                        $"{r.TotalCost / 1000:F1} млн руб, скважин с нарушениями {r.Wells.Count(w => w.Violations.Count > 0)} ({sw.Elapsed.TotalSeconds:F1} с)");
                    return r;
                }

                // 2) перебор верхней границы размера куста K при нижней границе MinWellsPerPad
                for (int K = Math.Max(p.StartMaxWellsPerPad, p.MinWellsPerPad); K <= p.EndMaxWellsPerPad; K += Math.Max(1, p.VariantStep))
                {
                    string name = $"{scName}; {p.MinWellsPerPad}–{K} скв./куст";
                    var (pads, est) = o.Optimize(K);
                    if (pads == null)
                    {
                        infeasible.Add(new InfeasibleVariant(name, scName, K, $"скважины нельзя разбить на кусты по {p.MinWellsPerPad}–{K} скв."));
                        Console.WriteLine($"{name}: нет решения");
                        continue;
                    }
                    AddVariant(name, scName, K, pads, est);
                }

                // 3) варианты с заданным общим числом кустов
                if (p.RunPadCountSweep)
                {
                    var (lo, hi) = o.PadCountRange(p.EndMaxWellsPerPad);
                    for (int np = lo; np <= hi; np++)
                    {
                        string name = $"{scName}; {np} кустов по {p.MinWellsPerPad}–{p.EndMaxWellsPerPad} скв.";
                        var (pads, est) = o.Optimize(p.EndMaxWellsPerPad, np);
                        if (pads == null) { infeasible.Add(new InfeasibleVariant(name, scName, p.EndMaxWellsPerPad, "нет разбиения на такое число кустов")); continue; }
                        AddVariant(name, scName + " (число кустов)", p.EndMaxWellsPerPad, pads, est);
                    }
                }
            }

            // рекомендуемый вариант: сценарий с исходными параметрами, минимум скважин с нарушениями, затем затрат
            var best = results.Where(r => r.Kind == "opt" && r.SameAs == null && r.Scenario == scenarios[0].name)
                              .OrderBy(r => r.Violations).ThenBy(r => r.TotalCost).FirstOrDefault()
                       ?? results.OrderBy(r => r.Violations).ThenBy(r => r.TotalCost).First();

            Reporting.PrintSummary(results, best);
            foreach (var x in infeasible) Console.WriteLine($"{x.Name}: нет решения — {x.Reason}");
            Reporting.WriteCsv(outDir, results, p, infeasible);
            string json = Reporting.BuildJson(results, best, data, map, p, Path.GetFileName(input), Path.GetFileName(mapPath), infeasible);
            File.WriteAllText(Path.Combine(outDir, "results.json"), json, new UTF8Encoding(false));
            Reporting.WriteHtml(Path.Combine(outDir, "report.html"), json);
            Console.WriteLine($"\nГотово за {sw.Elapsed.TotalSeconds:F1} с. Результаты: {Path.GetFullPath(outDir)} (report.html, summary.csv, pads.csv, wells.csv, trajectories/*.csv)");
            return 0;
        }
        catch (Exception ex)
        {
            Console.Error.WriteLine("Ошибка: " + ex.Message);
            Console.Error.WriteLine(ex.StackTrace);
            return 1;
        }
    }

    private static List<((double n, double e) center, List<int> wells)> ToPads(Optimizer o, List<PadAssignment> pads) =>
        pads.Select(x => ((o.SiteN[x.Site], o.SiteE[x.Site]), x.Wells.OrderBy(w => w).ToList())).ToList();

    private static bool SameSig(List<(int, int[])> a, List<(int, int[])> b)
    {
        if (a.Count != b.Count) return false;
        for (int i = 0; i < a.Count; i++)
            if (a[i].Item1 != b[i].Item1 || !a[i].Item2.SequenceEqual(b[i].Item2)) return false;
        return true;
    }


    private static string Arg(string[] a, string key)
    {
        for (int i = 0; i < a.Length - 1; i++) if (a[i] == key) return a[i + 1];
        return null;
    }
}
