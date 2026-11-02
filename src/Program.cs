using System.Diagnostics;
using System.Text;

namespace WellPadOptimizer;

/// <summary>
/// Кустование скважин и расчёт траекторий.
/// Запуск:  WellPadOptimizer --input data/input_mass.json [--map data/terrain_map.asc] [--settings settings.json] [--out results]
/// </summary>
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
            Optimizer lastOpt = null; Evaluator lastEval = null; Parameters lastPrm = null; string lastScenario = null;
            foreach (var (scName, prm) in scenarios)
            {
                var o = new Optimizer(prm, data, map) { SiteN = opt.SiteN, SiteE = opt.SiteE, Sand = opt.Sand };
                o.BuildCostMatrix();
                Console.WriteLine($"[{scName}] допустимых пар площадка–скважина: {CountFinite(o.C)}; отказы: " +
                    string.Join(", ", o.RejectStats.Where(k => k.Key != CheckResult.Ok).Select(k => $"{k.Key}={k.Value}")));
                var ev = new Evaluator(prm, data, map);
                List<(int, int[])> prevSig = null; VariantResult prev = null;
                for (int K = p.StartMaxWellsPerPad; K <= p.EndMaxWellsPerPad; K += Math.Max(1, p.VariantStep))
                {
                    string name = $"{scName}; ≤{K} скв./куст";
                    var (pads, est) = o.Optimize(K);
                    var sig = pads.Select(x => (x.Site, x.Wells.OrderBy(w => w).ToArray())).OrderBy(x => x.Site).ToList();
                    if (prev != null && SameSig(sig, prevSig))
                    {
                        results.Add(prev.AliasOf(name, K));
                        continue;
                    }
                    var r = ev.Evaluate(name, "opt", scName, K, ToPads(o, pads));
                    r.OptimizerEstimate = est;
                    results.Add(r); prev = r; prevSig = sig;
                    Console.WriteLine($"{name}: кустов {r.PadCount} [{string.Join(",", r.Pads.Select(x => x.WellCount))}], " +
                        $"{r.TotalCost / 1000:F1} млн руб, нарушений {r.Violations} ({sw.Elapsed.TotalSeconds:F1} с)");
                }
                lastOpt = o; lastEval = ev; lastPrm = prm; lastScenario = scName;
            }

            // 3) варианты с ограничением числа кустов (кривая «число кустов — затраты»)
            if (p.RunPadCountSweep && lastOpt != null)
            {
                int start = results.Where(r => r.Scenario == lastScenario && r.SameAs == null).Min(r => r.PadCount) - 1;
                for (int mp = start; mp >= 1; mp--)
                {
                    var (pads, est) = lastOpt.Optimize(p.EndMaxWellsPerPad, mp);
                    if (pads == null) { Console.WriteLine($"Не более {mp} кустов: решения нет"); break; }
                    string name = $"{lastScenario}; не более {mp} кустов";
                    var r = lastEval.Evaluate(name, "opt", "Лимит кустов", p.EndMaxWellsPerPad, ToPads(lastOpt, pads));
                    r.OptimizerEstimate = est;
                    results.Add(r);
                    Console.WriteLine($"{name}: [{string.Join(",", r.Pads.Select(x => x.WellCount))}], {r.TotalCost / 1000:F1} млн руб");
                }
            }

            // лучший допустимый вариант при исходных параметрах (сценарий 1)
            var best = results.Where(r => r.Kind == "opt" && r.Violations == 0 && r.SameAs == null && r.Scenario == scenarios[0].name)
                              .OrderBy(r => r.TotalCost).FirstOrDefault()
                       ?? results.OrderBy(r => r.Violations).ThenBy(r => r.TotalCost).First();

            Reporting.PrintSummary(results, best);
            Reporting.WriteCsv(outDir, results, p);
            string json = Reporting.BuildJson(results, best, data, map, p, Path.GetFileName(input), Path.GetFileName(mapPath));
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

    private static int CountFinite(double[,] c)
    {
        int n = 0;
        foreach (var v in c) if (!double.IsInfinity(v)) n++;
        return n;
    }

    private static string Arg(string[] a, string key)
    {
        for (int i = 0; i < a.Length - 1; i++) if (a[i] == key) return a[i + 1];
        return null;
    }
}
