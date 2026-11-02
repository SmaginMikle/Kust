using System.Globalization;
using System.Text;

namespace WellPadOptimizer;

/// <summary>
/// Карта желательности местности — растр ESRI ASCII Grid (.asc).
/// Коды классов задаются в настройках (по умолчанию 1 — суходол, 2 — болото, 3 — водоём, 9 — запрещено).
/// NODATA и точки вне карты считаются запрещёнными.
/// </summary>
public sealed class TerrainMap
{
    public int NCols { get; private set; }
    public int NRows { get; private set; }
    public double XLL { get; private set; }     // восток левого нижнего угла
    public double YLL { get; private set; }     // север левого нижнего угла
    public double CellSize { get; private set; }
    public int NoData { get; private set; } = -9999;
    public int ForbiddenCode { get; set; } = 9;
    private short[] _grid;                       // построчно, строка 0 — северная

    public int this[int row, int col] => _grid[row * NCols + col];

    /// <summary>Класс местности в точке (E — восток, N — север).</summary>
    public int ClassAt(double e, double n)
    {
        int ix = (int)Math.Floor((e - XLL) / CellSize);
        int iy = NRows - 1 - (int)Math.Floor((n - YLL) / CellSize);
        if (ix < 0 || ix >= NCols || iy < 0 || iy >= NRows) return ForbiddenCode;
        int c = _grid[iy * NCols + ix];
        return c == NoData ? ForbiddenCode : c;
    }

    public static TerrainMap LoadAsc(string path)
    {
        var map = new TerrainMap();
        using var sr = new StreamReader(path);
        var hdr = new Dictionary<string, double>(StringComparer.OrdinalIgnoreCase);
        string line;
        while (hdr.Count < 6 && (line = sr.ReadLine()) != null)
        {
            var p = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
            if (p.Length >= 2) hdr[p[0]] = double.Parse(p[1], CultureInfo.InvariantCulture);
        }
        map.NCols = (int)hdr["ncols"];
        map.NRows = (int)hdr["nrows"];
        map.XLL = hdr.ContainsKey("xllcorner") ? hdr["xllcorner"] : hdr["xllcenter"] - hdr["cellsize"] / 2;
        map.YLL = hdr.ContainsKey("yllcorner") ? hdr["yllcorner"] : hdr["yllcenter"] - hdr["cellsize"] / 2;
        map.CellSize = hdr["cellsize"];
        if (hdr.TryGetValue("nodata_value", out var nd)) map.NoData = (int)nd;
        map._grid = new short[map.NCols * map.NRows];
        int k = 0;
        while ((line = sr.ReadLine()) != null && k < map._grid.Length)
        {
            foreach (var tok in line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries))
            {
                if (k >= map._grid.Length) break;
                map._grid[k++] = (short)Math.Round(double.Parse(tok, CultureInfo.InvariantCulture));
            }
        }
        if (k != map._grid.Length) throw new InvalidDataException($"Карта {path}: прочитано {k} из {map._grid.Length} ячеек");
        return map;
    }

    public void SaveAsc(string path)
    {
        var sb = new StringBuilder();
        sb.AppendLine($"ncols {NCols}");
        sb.AppendLine($"nrows {NRows}");
        sb.AppendLine(FormattableString.Invariant($"xllcorner {XLL:F1}"));
        sb.AppendLine(FormattableString.Invariant($"yllcorner {YLL:F1}"));
        sb.AppendLine(FormattableString.Invariant($"cellsize {CellSize}"));
        sb.AppendLine($"NODATA_value {NoData}");
        for (int r = 0; r < NRows; r++)
        {
            for (int c = 0; c < NCols; c++) { if (c > 0) sb.Append(' '); sb.Append(_grid[r * NCols + c]); }
            sb.AppendLine();
        }
        File.WriteAllText(path, sb.ToString());
    }

    /// <summary>
    /// Синтетическая карта (если файл карты не задан): фрактальный шум → суходол/болото/водоём,
    /// плюс круговые запретные зоны и линейные коридоры. Детерминирована (seed).
    /// </summary>
    public static TerrainMap Synthetic(IReadOnlyList<WellTarget> wells, double cell = 50, double margin = 2500, int seed = 20261102)
    {
        double minE = wells.Min(w => Math.Min(w.T1.E, w.T3.E)) - margin, maxE = wells.Max(w => Math.Max(w.T1.E, w.T3.E)) + margin;
        double minN = wells.Min(w => Math.Min(w.T1.N, w.T3.N)) - margin, maxN = wells.Max(w => Math.Max(w.T1.N, w.T3.N)) + margin;
        double x0 = Math.Floor(minE / 1000) * 1000, y0 = Math.Floor(minN / 1000) * 1000;
        double x1 = Math.Ceiling(maxE / 1000) * 1000, y1 = Math.Ceiling(maxN / 1000) * 1000;
        var map = new TerrainMap
        {
            NCols = (int)((x1 - x0) / cell), NRows = (int)((y1 - y0) / cell),
            XLL = x0, YLL = y0, CellSize = cell, ForbiddenCode = 9,
        };
        int nc = map.NCols, nr = map.NRows;
        var rnd = new Random(seed);
        var f = new double[nr * nc];
        void Octave(double scale, double amp)
        {
            int gy = (int)(nr * cell / scale) + 3, gx = (int)(nc * cell / scale) + 3;
            var g = new double[gy, gx];
            for (int i = 0; i < gy; i++) for (int j = 0; j < gx; j++) g[i, j] = rnd.NextDouble();
            for (int r = 0; r < nr; r++)
            {
                double yy = (r + 0.5) * cell / scale; int iy = (int)yy; double fy = yy - iy; fy = fy * fy * (3 - 2 * fy);
                for (int c = 0; c < nc; c++)
                {
                    double xx = (c + 0.5) * cell / scale; int ix = (int)xx; double fx = xx - ix; fx = fx * fx * (3 - 2 * fx);
                    double v = (g[iy, ix] * (1 - fx) + g[iy, ix + 1] * fx) * (1 - fy) + (g[iy + 1, ix] * (1 - fx) + g[iy + 1, ix + 1] * fx) * fy;
                    f[r * nc + c] += amp * v;
                }
            }
        }
        Octave(3000, 0.5); Octave(1200, 0.3); Octave(500, 0.15); Octave(200, 0.05);
        var sorted = (double[])f.Clone(); Array.Sort(sorted);
        double q50 = sorted[(int)(0.50 * (sorted.Length - 1))], q88 = sorted[(int)(0.88 * (sorted.Length - 1))];
        map._grid = new short[nr * nc];
        for (int i = 0; i < f.Length; i++) map._grid[i] = (short)(f[i] > q88 ? 3 : f[i] > q50 ? 2 : 1);

        // запретные зоны: окружности возле целей
        for (int k = 0; k < 14; k++)
        {
            var w = wells[rnd.Next(wells.Count)];
            double ce = w.T1.E + Gauss(rnd) * 900, cn = w.T1.N + Gauss(rnd) * 900, rad = 250 + rnd.NextDouble() * 300;
            for (int r = 0; r < nr; r++)
                for (int c = 0; c < nc; c++)
                {
                    double e = x0 + (c + 0.5) * cell, n = y1 - (r + 0.5) * cell;
                    if ((e - ce) * (e - ce) + (n - cn) * (n - cn) < rad * rad) map._grid[r * nc + c] = 9;
                }
        }
        return map;
    }

    private static double Gauss(Random r)
    {
        double u1 = 1.0 - r.NextDouble(), u2 = r.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
    }

    /// <summary>Вырезка растра (для отчёта) с прореживанием step.</summary>
    public int[][] Crop(double e0, double n0, double e1, double n1, out double cropE0, out double cropN1, out double outCell, int step = 1)
    {
        int c0 = Math.Max(0, (int)Math.Floor((e0 - XLL) / CellSize));
        int c1 = Math.Min(NCols - 1, (int)Math.Floor((e1 - XLL) / CellSize));
        int r0 = Math.Max(0, NRows - 1 - (int)Math.Floor((n1 - YLL) / CellSize));
        int r1 = Math.Min(NRows - 1, NRows - 1 - (int)Math.Floor((n0 - YLL) / CellSize));
        cropE0 = XLL + c0 * CellSize;
        cropN1 = YLL + (NRows - r0) * CellSize;
        outCell = CellSize * step;
        var rows = new List<int[]>();
        for (int r = r0; r <= r1; r += step)
        {
            var row = new List<int>();
            for (int c = c0; c <= c1; c += step)
            {
                int v = _grid[r * NCols + c];
                row.Add(v == NoData ? ForbiddenCode : v);
            }
            rows.Add(row.ToArray());
        }
        return rows.ToArray();
    }
}
