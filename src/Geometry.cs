namespace WellPadOptimizer;

/// <summary>Точка/вектор в системе (N — север, E — восток, Z — вертикальная глубина вниз).</summary>
public readonly struct Vec3
{
    public readonly double N, E, Z;
    public Vec3(double n, double e, double z) { N = n; E = e; Z = z; }

    public static readonly Vec3 Down = new(0, 0, 1);
    public static Vec3 operator +(Vec3 a, Vec3 b) => new(a.N + b.N, a.E + b.E, a.Z + b.Z);
    public static Vec3 operator -(Vec3 a, Vec3 b) => new(a.N - b.N, a.E - b.E, a.Z - b.Z);
    public static Vec3 operator *(double k, Vec3 a) => new(k * a.N, k * a.E, k * a.Z);
    public static Vec3 operator *(Vec3 a, double k) => new(k * a.N, k * a.E, k * a.Z);
    public double Dot(Vec3 b) => N * b.N + E * b.E + Z * b.Z;
    public double Length => Math.Sqrt(N * N + E * E + Z * Z);
    public double HorizontalLength => Math.Sqrt(N * N + E * E);
    public Vec3 Unit() { double l = Length; return l < 1e-15 ? this : new Vec3(N / l, E / l, Z / l); }
    public override string ToString() => $"({N:F2}; {E:F2}; {Z:F2})";
}

public static class Geo
{
    public const double Deg = Math.PI / 180.0;

    /// <summary>Азимут (°, от севера по часовой) по приращениям dN, dE.</summary>
    public static double Azimuth(double dN, double dE)
    {
        double a = Math.Atan2(dE, dN) / Deg;
        return (a + 360.0) % 360.0;
    }

    public static double AngleDiff(double a, double b)
    {
        double d = Math.Abs(a - b) % 360.0;
        return Math.Min(d, 360.0 - d);
    }

    public static double Clamp1(double x) => x < -1 ? -1 : (x > 1 ? 1 : x);
    public static double Hypot(double a, double b) => Math.Sqrt(a * a + b * b);
}
