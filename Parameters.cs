using System.Text.Json;

namespace WellPadOptimizer;

/// <summary>Группа стоимости строительства куста (тыс. руб. на 1 скважину).</summary>
public sealed class PadCostGroup
{
    public int MaxWells { get; set; }
    public double CostPerWell { get; set; }
}

/// <summary>Класс местности на карте желательности.</summary>
public sealed class TerrainClass
{
    public int Code { get; set; }
    public string Name { get; set; } = "";
    /// <summary>Высота отсыпки, м.</summary>
    public double FillHeight { get; set; }
}

/// <summary>
/// Параметры расчёта. Значения по умолчанию соответствуют окну настроек
/// (максимальный отход 1500 м, MD 4000 м, темпы набора 0,5 и 1,5 °/10 м и т.д.).
/// </summary>
public sealed class Parameters
{
    // --- пользовательские параметры ---
    public double MaxOffset { get; set; } = 1500;          // Максимальный отход от устья, м
    public double MaxMD { get; set; } = 4000;              // Максимальная длина по стволу, м
    public bool LimitsToBottomHole { get; set; } = false;  // «На забой»: ограничения до T3, иначе до T1
    public double Dls1 { get; set; } = 0.5;                // Темп первого набора, °/DlsBaseLength
    public double Dls2 { get; set; } = 1.5;                // Темп второго набора, °/DlsBaseLength
    public double DlsBaseLength { get; set; } = 10;        // База темпа, м
    public bool UseAzimuthTurn { get; set; } = true;       // Учитывать максимальный разворот
    public double MaxAzimuthTurn { get; set; } = 90;       // Максимальный разворот по азимуту, °
    public double VerticalLength { get; set; } = 40;       // Длина вертикального участка, м
    public int VariantStep { get; set; } = 1;              // Пропуск вариантов
    public int StartMaxWellsPerPad { get; set; } = 10;     // Начальное макс. число скважин на кусте
    public int EndMaxWellsPerPad { get; set; } = 24;       // Конечное макс. число скважин на кусте
    public int MinWellsPerPad { get; set; } = 10;          // мин. число скважин на кусте (жёсткое ограничение)
    // Поиск площадок до SoftRadiusFactor·MaxOffset: технологические ограничения могут нарушаться со штрафом,
    // если при заданном размере куста иначе решения нет. Штрафы в тыс. руб.
    public double SoftRadiusFactor { get; set; } = 1.6;
    public double PenaltyViolation { get; set; } = 1.0e6;
    public double PenaltyOffsetPerM { get; set; } = 1000;
    public double PenaltyAzimuthPerDeg { get; set; } = 10000;
    public double PenaltyMdPerM { get; set; } = 1000;
    public double MaxHookLoad { get; set; } = 180;         // тс
    public double MaxTorque { get; set; } = 54;            // кН·м

    // --- стоимостные параметры (тыс. руб.) ---
    public double CostNNS { get; set; } = 65;              // за метр ННС (до T1)
    public double CostGS { get; set; } = 55;               // за метр ГС (T1–T3)
    public List<PadCostGroup> PadCostGroups { get; set; } = new()
    {
        new PadCostGroup { MaxWells = 12, CostPerWell = 16000 },
        new PadCostGroup { MaxWells = 19, CostPerWell = 13000 },
        new PadCostGroup { MaxWells = 24, CostPerWell = 12000 },
    };
    public double SandCostPer1000m3 { get; set; } = 1550;

    // --- расстановка устьев ---
    public int WellsInGroup { get; set; } = 4;
    public double DistanceInGroup { get; set; } = 6;
    public double DistanceBetweenGroups { get; set; } = 15;

    // --- геометрия и отсыпка площадки ---
    public double PadWidth { get; set; } = 120;
    public double PadEndMargin { get; set; } = 60;
    public double PadMinLength { get; set; } = 100;
    public double EmbankmentSlope { get; set; } = 1.5;     // заложение откоса 1:m
    public List<TerrainClass> Terrain { get; set; } = new()
    {
        new TerrainClass { Code = 1, Name = "Суходол", FillHeight = 1.0 },
        new TerrainClass { Code = 2, Name = "Болото", FillHeight = 2.5 },
        new TerrainClass { Code = 3, Name = "Водоём", FillHeight = 4.0 },
    };
    public int ForbiddenCode { get; set; } = 9;
    public double ForbiddenBuffer { get; set; } = 0;       // Мин. расстояние до запрещённой зоны, м
    public double ExtraPadGap { get; set; } = 0;           // Доп. расстояние между кустами, м

    // --- дискретизация поиска ---
    public double SiteGridStep { get; set; } = 100;
    public double FootprintSampleStep { get; set; } = 10;
    public List<double> PadOrientations { get; set; } = new() { 0, 30, 60, 90, 120, 150 };

    // --- модель «мягкой» колонны для нагрузки на крюке и момента ---
    public double PipeWeightAir { get; set; } = 0.284;     // кН/м (127 мм БТ)
    public double MudDensity { get; set; } = 1.20;         // г/см3
    public double SteelDensity { get; set; } = 7.85;
    public double FrictionCased { get; set; } = 0.25;
    public double FrictionOpenHole { get; set; } = 0.30;
    public double PipeRadius { get; set; } = 0.08;         // м (по замкам)

    // --- траектории и антиколлизия ---
    public double TrajectoryStep { get; set; } = 10;
    public double MinSeparation { get; set; } = 8;         // цель при раскладке скважин по слотам
    public double AntiCollisionMaxMD { get; set; } = 1000;
    public double SeparationWarning { get; set; } = 5;

    // --- набор вариантов ---
    public bool RunScenarioWithoutAzimuthTurn { get; set; } = true;
    public bool RunPadCountSweep { get; set; } = true;
    public string PadPrefix { get; set; } = "К-";
    public int StartPadNumber { get; set; } = 1;

    public double R1 => DlsBaseLength / (Dls1 * Math.PI / 180.0);
    public double R2 => DlsBaseLength / (Dls2 * Math.PI / 180.0);
    public int MaxPadCapacity => PadCostGroups.Max(g => g.MaxWells);

    public Parameters Clone() => (Parameters)MemberwiseClone();

    public static Parameters Load(string path)
    {
        if (string.IsNullOrEmpty(path) || !File.Exists(path)) return new Parameters();
        var opt = new JsonSerializerOptions
        {
            PropertyNameCaseInsensitive = true,
            ReadCommentHandling = JsonCommentHandling.Skip,
            AllowTrailingCommas = true,
        };
        return JsonSerializer.Deserialize<Parameters>(File.ReadAllText(path), opt) ?? new Parameters();
    }

    // ---------------- функции стоимости и геометрии площадки ----------------

    /// <summary>Стоимость строительства куста (тыс. руб.) в зависимости от числа скважин.</summary>
    public double PadBuildCost(int n)
    {
        if (n <= 0) return 0;
        foreach (var g in PadCostGroups.OrderBy(g => g.MaxWells))
            if (n <= g.MaxWells) return g.CostPerWell * n;
        return double.PositiveInfinity;
    }

    /// <summary>Длина ряда устьев: группы по WellsInGroup с шагом DistanceInGroup, между группами DistanceBetweenGroups.</summary>
    public double RowLength(int n)
    {
        if (n <= 1) return 0;
        int groups = (n - 1) / WellsInGroup;
        return (n - 1) * DistanceInGroup + groups * (DistanceBetweenGroups - DistanceInGroup);
    }

    /// <summary>Смещения слотов вдоль оси ряда относительно центра куста.</summary>
    public double[] SlotOffsets(int n)
    {
        var s = new double[n];
        for (int k = 0; k < n; k++)
            s[k] = k * DistanceInGroup + (k / WellsInGroup) * (DistanceBetweenGroups - DistanceInGroup);
        double c = n > 0 ? s[n - 1] / 2 : 0;
        for (int k = 0; k < n; k++) s[k] -= c;
        return s;
    }

    public double PadLength(int n) => Math.Max(PadMinLength, RowLength(n) + 2 * PadEndMargin);
    public double PadDiagonal(int n) => Math.Sqrt(PadLength(n) * PadLength(n) + PadWidth * PadWidth);

    /// <summary>Объём отсыпки (м3) по долям классов местности под площадкой (призматоид с откосами).</summary>
    public double SandVolume(IReadOnlyDictionary<int, double> fractions, int n)
    {
        double L = PadLength(n), W = PadWidth, m = EmbankmentSlope, v = 0;
        foreach (var t in Terrain)
        {
            fractions.TryGetValue(t.Code, out double f);
            double h = t.FillHeight;
            v += f * h * (L + m * h) * (W + m * h);
        }
        return v;
    }

    public double SandCost(double volume) => volume / 1000.0 * SandCostPer1000m3;
}
