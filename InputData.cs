using System.Text.Json;

namespace WellPadOptimizer;

/// <summary>Цель скважины: T1 — начало горизонтального участка, T3 — окончание.</summary>
public sealed class WellTarget
{
    public int Index { get; init; }
    public string Name { get; init; } = "";
    public string Uid { get; init; } = "";
    public Vec3 T1 { get; init; }
    public Vec3 T3 { get; init; }
    public double HorizontalLength => (T3 - T1).Length;
    public bool HasHorizontal => HorizontalLength > 1e-6;
    public string Key => $"{T1.N:F2}|{T1.E:F2}|{T3.N:F2}|{T3.E:F2}";
}

/// <summary>Фактическая (существующая) скважина из trajRealWells.</summary>
public sealed class RealWell
{
    public string Name { get; init; } = "";
    public double X { get; init; }   // восток
    public double Y { get; init; }   // север
    public double Altitude { get; init; }
    public List<Vec3> Points { get; } = new();
}

/// <summary>Вариант кустования, заданный во входном файле (устья заданы).</summary>
public sealed class InputVariant
{
    public int Id { get; init; }
    /// <summary>pad_id → список (индекс цели, устье N, устье E).</summary>
    public SortedDictionary<int, List<(int well, double n, double e)>> Pads { get; } = new();
}

public sealed class InputData
{
    public List<WellTarget> Wells { get; } = new();
    public List<RealWell> RealWells { get; } = new();
    public List<InputVariant> Variants { get; } = new();
    public int ProjectedCrs { get; private set; }

    /// <summary>
    /// Чтение формата rndrilling-v3: massiveDataVariants[].trajBuildRequest[] с T1/T3 = [North, East, TVD].
    /// Набор целей берётся из первого варианта; остальные варианты сопоставляются по координатам T1/T3.
    /// </summary>
    public static InputData Load(string path)
    {
        var data = new InputData();
        using var doc = JsonDocument.Parse(File.ReadAllBytes(path),
            new JsonDocumentOptions { AllowTrailingCommas = true });
        var root = doc.RootElement;

        if (root.TryGetProperty("mapDefinition", out var md) && md.TryGetProperty("projected_crs", out var crs))
            data.ProjectedCrs = crs.GetInt32();

        var variants = root.GetProperty("massiveDataVariants");
        var keyToIndex = new Dictionary<string, int>();
        int vi = 0;
        foreach (var v in variants.EnumerateArray())
        {
            var iv = new InputVariant { Id = vi };
            foreach (var w in v.GetProperty("trajBuildRequest").EnumerateArray())
            {
                var t1 = ReadPoint(w.GetProperty("T1"));
                var t3 = ReadPoint(w.GetProperty("T3"));
                string key = $"{t1.N:F2}|{t1.E:F2}|{t3.N:F2}|{t3.E:F2}";
                if (!keyToIndex.TryGetValue(key, out int idx))
                {
                    idx = data.Wells.Count;
                    keyToIndex[key] = idx;
                    data.Wells.Add(new WellTarget
                    {
                        Index = idx,
                        Name = $"W{idx + 1:00}",
                        Uid = w.TryGetProperty("uid", out var u) ? u.GetString() : idx.ToString(),
                        T1 = t1,
                        T3 = t3,
                    });
                }
                int pad = w.TryGetProperty("pad_id", out var p) ? p.GetInt32() : 0;
                if (w.TryGetProperty("wh_north", out var wn) && w.TryGetProperty("wh_east", out var we)
                    && wn.ValueKind == JsonValueKind.Number && we.ValueKind == JsonValueKind.Number)
                {
                    if (!iv.Pads.TryGetValue(pad, out var list)) iv.Pads[pad] = list = new();
                    list.Add((idx, wn.GetDouble(), we.GetDouble()));
                }
            }
            if (iv.Pads.Count > 0) data.Variants.Add(iv);
            vi++;
        }

        if (root.TryGetProperty("trajRealWells", out var rw))
        {
            foreach (var w in rw.EnumerateArray())
            {
                var r = new RealWell
                {
                    Name = w.GetProperty("Name").GetString(),
                    X = w.GetProperty("CoordsX").GetDouble(),
                    Y = w.GetProperty("CoordsY").GetDouble(),
                    Altitude = w.TryGetProperty("Altitude", out var a) ? a.GetDouble() : 0,
                };
                var xs = w.GetProperty("X"); var ys = w.GetProperty("Y"); var zs = w.GetProperty("AbsoluteDepth");
                int n = Math.Min(xs.GetArrayLength(), Math.Min(ys.GetArrayLength(), zs.GetArrayLength()));
                for (int i = 0; i < n; i++)
                    r.Points.Add(new Vec3(r.Y + ys[i].GetDouble(), r.X + xs[i].GetDouble(), zs[i].GetDouble()));
                data.RealWells.Add(r);
            }
        }
        return data;
    }

    private static Vec3 ReadPoint(JsonElement a) => new(a[0].GetDouble(), a[1].GetDouble(), a[2].GetDouble());
}
