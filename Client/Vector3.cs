namespace GlxGears.Client;

// WebGL Programming Guide: Interactive 3D Graphics Programming with WebGL 
//
// https://sites.google.com/site/webglbook/

public class Vector3
{
    public float[] Elements { get; } = new float[3];

    public Vector3(params float[] victor)
    {
        foreach (int idx in Enumerable.Range(0, victor.Length))
        {
            Elements[idx] = victor[idx];
        }
    }

    public Vector3 Normalize()
    {
        float[] v = Elements;

        float c = v[0], d = v[1], e = v[2];

        float g = (float)Math.Sqrt(c*c+d*d+e*e);

        if (g == 1) {
            return this;
        }

        g = 1/g;

        v[0] = c*g; v[1] = d*g; v[2] = e*g;

        return this;
    }
}
