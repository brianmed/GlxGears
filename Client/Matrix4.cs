namespace GlxGears.Client;

// WebGL Programming Guide: Interactive 3D Graphics Programming with WebGL 
//
// https://sites.google.com/site/webglbook/

public class Matrix4
{
    public float[] Elements { get; } = new float[16];

    public Matrix4()
    {
        SetIdentity();
    }

    public Matrix4 Invert()
    {
        return this.SetInverseOf(this);
    }

    public Matrix4 Multiply(Matrix4 other)
    {
        float ai0, ai1, ai2, ai3;
        
        // Calculate e = a * b
        float[] e = Elements;
        float[] a = Elements;
        float[] b = other.Elements;
        
        // If e equals b, copy b to temporary matrix.
        if (e == b) {
            b = new float[16];
            
            for (int i = 0; i < 16; ++i)
            {
                b[i] = e[i];
            }
        }
        
        for (int i = 0; i < 4; i++)
        {
            ai0=a[i];  ai1=a[i+4];  ai2=a[i+8];  ai3=a[i+12];
            e[i]    = ai0 * b[0]  + ai1 * b[1]  + ai2 * b[2]  + ai3 * b[3];
            e[i+4]  = ai0 * b[4]  + ai1 * b[5]  + ai2 * b[6]  + ai3 * b[7];
            e[i+8]  = ai0 * b[8]  + ai1 * b[9]  + ai2 * b[10] + ai3 * b[11];
            e[i+12] = ai0 * b[12] + ai1 * b[13] + ai2 * b[14] + ai3 * b[15];
        }
        
        return this;
    }

    public Matrix4 Rotate(float angle, float x, float y, float z)
    {
        Matrix4 rotate = new Matrix4()
            .SetRotate(angle, x, y, z);

        return Multiply(rotate);
    }

    /**
     * Multiply the matrix for scaling from the right.
     * @param x The scale factor along the X axis
     * @param y The scale factor along the Y axis
     * @param z The scale factor along the Z axis
     * @return this
     */
    public Matrix4 Scale(float x, float y, float z)
    {
        float[] e = Elements;
        
        e[0] *= x;  e[4] *= y;  e[8]  *= z;
        e[1] *= x;  e[5] *= y;  e[9]  *= z;
        e[2] *= x;  e[6] *= y;  e[10] *= z;
        e[3] *= x;  e[7] *= y;  e[11] *= z;
        
        return this;
    }

    /**
     * Set the matrix for scaling.
     * @param x The scale factor along the X axis
     * @param y The scale factor along the Y axis
     * @param z The scale factor along the Z axis
     * @return this
     */
    public Matrix4 SetScale(float x, float y, float z)
    {
        float[] e = Elements;
        
        e[0] = x;  e[4] = 0;  e[8]  = 0;  e[12] = 0;
        e[1] = 0;  e[5] = y;  e[9]  = 0;  e[13] = 0;
        e[2] = 0;  e[6] = 0;  e[10] = z;  e[14] = 0;
        e[3] = 0;  e[7] = 0;  e[11] = 0;  e[15] = 1;

        return this;
    }

    /**
     * Set the viewing matrix.
     * @param eyeX, eyeY, eyeZ The position of the eye point.
     * @param atX, atY, atZ The position of the reference point.
     * @param upX, upY, upZ The direction of the up vector.
     * @return this
     */
    public Matrix4 SetLookAt(float eyeX, float eyeY, float eyeZ, float atX, float atY, float atZ, float upX, float upY, float upZ)
    {
        // fx, fy, fz, rlf, sx, sy, sz, rls, ux, uy, uz;
        
        float fx = atX - eyeX;
        float fy = atY - eyeY;
        float fz = atZ - eyeZ;
        
        // Normalize f.
        float rlf = (float)(1.0 / Math.Sqrt(fx*fx + fy*fy + fz*fz));
        fx *= rlf;
        fy *= rlf;
        fz *= rlf;
        
        // Calculate cross product of f and up.
        float sx = fy * upZ - fz * upY;
        float sy = fz * upX - fx * upZ;
        float sz = fx * upY - fy * upX;
        
        // Normalize s.
        float rls = (float)(1.0 / Math.Sqrt(sx*sx + sy*sy + sz*sz));
        sx *= rls;
        sy *= rls;
        sz *= rls;
        
        // Calculate cross product of s and f.
        float ux = sy * fz - sz * fy;
        float uy = sz * fx - sx * fz;
        float uz = sx * fy - sy * fx;
        
        // Set to this.
        float[] e = Elements;
        e[0] = sx;
        e[1] = ux;
        e[2] = -fx;
        e[3] = 0;
        
        e[4] = sy;
        e[5] = uy;
        e[6] = -fy;
        e[7] = 0;
        
        e[8] = sz;
        e[9] = uz;
        e[10] = -fz;
        e[11] = 0;
        
        e[12] = 0;
        e[13] = 0;
        e[14] = 0;
        e[15] = 1;
        
        // Translate.
        return Translate(-eyeX, -eyeY, -eyeZ);
    }

    public Matrix4 SetIdentity()
    {
        float[] e = Elements;
        
        e[0] = 1;   e[4] = 0;   e[8]  = 0;   e[12] = 0;
        e[1] = 0;   e[5] = 1;   e[9]  = 0;   e[13] = 0;
        e[2] = 0;   e[6] = 0;   e[10] = 1;   e[14] = 0;
        e[3] = 0;   e[7] = 0;   e[11] = 0;   e[15] = 1;
        
        return this;
    }


    /**
     * Calculate the inverse matrix of specified matrix, and set to this.
     * @param other The source matrix
     * @return this
     */
    public Matrix4 SetInverseOf(Matrix4 other)
    {
        float[] s = other.Elements;

        float[] d = Elements;

        float[] inv = new float[16];
        
        inv[0]  =   s[5]*s[10]*s[15] - s[5] *s[11]*s[14] - s[9] *s[6]*s[15]
                  + s[9]*s[7] *s[14] + s[13]*s[6] *s[11] - s[13]*s[7]*s[10];
        inv[4]  = - s[4]*s[10]*s[15] + s[4] *s[11]*s[14] + s[8] *s[6]*s[15]
                  - s[8]*s[7] *s[14] - s[12]*s[6] *s[11] + s[12]*s[7]*s[10];
        inv[8]  =   s[4]*s[9] *s[15] - s[4] *s[11]*s[13] - s[8] *s[5]*s[15]
                  + s[8]*s[7] *s[13] + s[12]*s[5] *s[11] - s[12]*s[7]*s[9];
        inv[12] = - s[4]*s[9] *s[14] + s[4] *s[10]*s[13] + s[8] *s[5]*s[14]
                  - s[8]*s[6] *s[13] - s[12]*s[5] *s[10] + s[12]*s[6]*s[9];
        
        inv[1]  = - s[1]*s[10]*s[15] + s[1] *s[11]*s[14] + s[9] *s[2]*s[15]
                  - s[9]*s[3] *s[14] - s[13]*s[2] *s[11] + s[13]*s[3]*s[10];
        inv[5]  =   s[0]*s[10]*s[15] - s[0] *s[11]*s[14] - s[8] *s[2]*s[15]
                  + s[8]*s[3] *s[14] + s[12]*s[2] *s[11] - s[12]*s[3]*s[10];
        inv[9]  = - s[0]*s[9] *s[15] + s[0] *s[11]*s[13] + s[8] *s[1]*s[15]
                  - s[8]*s[3] *s[13] - s[12]*s[1] *s[11] + s[12]*s[3]*s[9];
        inv[13] =   s[0]*s[9] *s[14] - s[0] *s[10]*s[13] - s[8] *s[1]*s[14]
                  + s[8]*s[2] *s[13] + s[12]*s[1] *s[10] - s[12]*s[2]*s[9];
        
        inv[2]  =   s[1]*s[6]*s[15] - s[1] *s[7]*s[14] - s[5] *s[2]*s[15]
                  + s[5]*s[3]*s[14] + s[13]*s[2]*s[7]  - s[13]*s[3]*s[6];
        inv[6]  = - s[0]*s[6]*s[15] + s[0] *s[7]*s[14] + s[4] *s[2]*s[15]
                  - s[4]*s[3]*s[14] - s[12]*s[2]*s[7]  + s[12]*s[3]*s[6];
        inv[10] =   s[0]*s[5]*s[15] - s[0] *s[7]*s[13] - s[4] *s[1]*s[15]
                  + s[4]*s[3]*s[13] + s[12]*s[1]*s[7]  - s[12]*s[3]*s[5];
        inv[14] = - s[0]*s[5]*s[14] + s[0] *s[6]*s[13] + s[4] *s[1]*s[14]
                  - s[4]*s[2]*s[13] - s[12]*s[1]*s[6]  + s[12]*s[2]*s[5];
        
        inv[3]  = - s[1]*s[6]*s[11] + s[1]*s[7]*s[10] + s[5]*s[2]*s[11]
                  - s[5]*s[3]*s[10] - s[9]*s[2]*s[7]  + s[9]*s[3]*s[6];
        inv[7]  =   s[0]*s[6]*s[11] - s[0]*s[7]*s[10] - s[4]*s[2]*s[11]
                  + s[4]*s[3]*s[10] + s[8]*s[2]*s[7]  - s[8]*s[3]*s[6];
        inv[11] = - s[0]*s[5]*s[11] + s[0]*s[7]*s[9]  + s[4]*s[1]*s[11]
                  - s[4]*s[3]*s[9]  - s[8]*s[1]*s[7]  + s[8]*s[3]*s[5];
        inv[15] =   s[0]*s[5]*s[10] - s[0]*s[6]*s[9]  - s[4]*s[1]*s[10]
                  + s[4]*s[2]*s[9]  + s[8]*s[1]*s[6]  - s[8]*s[2]*s[5];
        
        float det = s[0]*inv[0] + s[1]*inv[4] + s[2]*inv[8] + s[3]*inv[12];

        if (det == 0) {
            return this;
        }

        det = 1 / det;

        for (int i = 0; i < 16; i++) {
            d[i] = inv[i] * det;
        }

        return this;
    }

    public Matrix4 SetFrustum(float left, float right, float bottom, float top, float near, float far)
    {
        float[] e = Elements;

        e[0]  = (2.0f * near) / (right - left);
        e[1]  = 0.0f;
        e[2]  = 0.0f;
        e[3]  = 0.0f;
        
        e[4]  = 0.0f;
        e[5]  = (2.0f * near) / (top - bottom);
        e[6]  = 0.0f;
        e[7]  = 0.0f;
        
        e[8]  = (right + left) / (right - left);
        e[9]  = (top + bottom) / (top - bottom);
        e[10] = -1.0f * ((far + near) / (far - near));
        e[11] = 0.0f;
        
        e[12] = 0.0f;
        e[13] = 0.0f;
        e[14] = -1.0f * ((2.0f * far * near) / far - near);
        e[15] = 0.0f;

        return this;
    }

    /**
     * Set the orthographic projection matrix.
     * @param left The coordinate of the left of clipping plane.
     * @param right The coordinate of the right of clipping plane.
     * @param bottom The coordinate of the bottom of clipping plane.
     * @param top The coordinate of the top top clipping plane.
     * @param near The distances to the nearer depth clipping plane. This value is minus if the plane is to be behind the viewer.
     * @param far The distances to the farther depth clipping plane. This value is minus if the plane is to be behind the viewer.
     * @return this
     */
    public Matrix4 SetOrtho(float left, float right, float bottom, float top, float near, float far)
    {
        float rw, rh, rd;
        
        rw = 1 / (right - left);
        rh = 1 / (top - bottom);
        rd = 1 / (far - near);
        
        float[] e = Elements;
        
        e[0]  = 2 * rw;
        e[1]  = 0;
        e[2]  = 0;
        e[3]  = 0;
        
        e[4]  = 0;
        e[5]  = 2 * rh;
        e[6]  = 0;
        e[7]  = 0;
        
        e[8]  = 0;
        e[9]  = 0;
        e[10] = -2 * rd;
        e[11] = 0;
        
        e[12] = -(right + left) * rw;
        e[13] = -(top + bottom) * rh;
        e[14] = -(far + near) * rd;
        e[15] = 1;
        
        return this;
    }

    /**
     * Set the perspective projection matrix by fovy and aspect.
     * @param fovy The angle between the upper and lower sides of the frustum.
     * @param aspect The aspect ratio of the frustum. (width/height)
     * @param near The distances to the nearer depth clipping plane. This value must be plus value.
     * @param far The distances to the farther depth clipping plane. This value must be plus value.
     * @return this
     */
    public Matrix4 SetPerspective(float fovy, float aspect, float near, float far)
    {
        float rd, s, ct;
        
        if (near == far || aspect == 0) {
          throw new Exception("null frustum");
        }
        if (near <= 0) {
          throw new Exception("near <= 0");
        }
        if (far <= 0) {
          throw new Exception("far <= 0");
        }
        
        fovy = (float)(Math.PI * fovy / 180.0 / 2.0);

        s = (float)Math.Sin(fovy);

        if (s == 0) {
          throw new Exception("null frustum");
        }
        
        rd = 1 / (far - near);
        ct = (float)(Math.Cos(fovy) / s);
        
        float[] e = Elements;
        
        e[0]  = ct / aspect;
        e[1]  = 0;
        e[2]  = 0;
        e[3]  = 0;
        
        e[4]  = 0;
        e[5]  = ct;
        e[6]  = 0;
        e[7]  = 0;
        
        e[8]  = 0;
        e[9]  = 0;
        e[10] = -(far + near) * rd;
        e[11] = -1;
        
        e[12] = 0;
        e[13] = 0;
        e[14] = -2 * near * far * rd;
        e[15] = 0;
        
        return this;
    }

    public Matrix4 SetRotate(float _angle, float x, float y, float z)
    {
        float angle = _angle;

        angle = (float)(Math.PI * angle / 180.0);

        float[] e = this.Elements;

        float s = (float)Math.Sin(angle);
        float c = (float)Math.Cos(angle);

        if (0 != x && 0 == y && 0 == z) {
          // Rotation around X axis
          if (x < 0) {
            s = -s;
          }
          e[0] = 1;  e[4] = 0;  e[ 8] = 0;  e[12] = 0;
          e[1] = 0;  e[5] = c;  e[ 9] =-s;  e[13] = 0;
          e[2] = 0;  e[6] = s;  e[10] = c;  e[14] = 0;
          e[3] = 0;  e[7] = 0;  e[11] = 0;  e[15] = 1;
        } else if (0 == x && 0 != y && 0 == z) {
          // Rotation around Y axis
          if (y < 0) {
            s = -s;
          }
          e[0] = c;  e[4] = 0;  e[ 8] = s;  e[12] = 0;
          e[1] = 0;  e[5] = 1;  e[ 9] = 0;  e[13] = 0;
          e[2] =-s;  e[6] = 0;  e[10] = c;  e[14] = 0;
          e[3] = 0;  e[7] = 0;  e[11] = 0;  e[15] = 1;
        } else if (0 == x && 0 == y && 0 != z) {
          // Rotation around Z axis
          if (z < 0) {
            s = -s;
          }
          e[0] = c;  e[4] =-s;  e[ 8] = 0;  e[12] = 0;
          e[1] = s;  e[5] = c;  e[ 9] = 0;  e[13] = 0;
          e[2] = 0;  e[6] = 0;  e[10] = 1;  e[14] = 0;
          e[3] = 0;  e[7] = 0;  e[11] = 0;  e[15] = 1;
        } else {
          // Rotation around another axis
          float len = (float)Math.Sqrt(x*x + y*y + z*z);
          if (len != 1) {
            float rlen = 1 / len;
            x *= rlen;
            y *= rlen;
            z *= rlen;
          }
          float nc = 1 - c;
          float xy = x * y;
          float yz = y * z;
          float zx = z * x;
          float xs = x * s;
          float ys = y * s;
          float zs = z * s;

          e[ 0] = x*x*nc +  c;
          e[ 1] = xy *nc + zs;
          e[ 2] = zx *nc - ys;
          e[ 3] = 0;

          e[ 4] = xy *nc - zs;
          e[ 5] = y*y*nc +  c;
          e[ 6] = yz *nc + xs;
          e[ 7] = 0;

          e[ 8] = zx *nc + ys;
          e[ 9] = yz *nc - xs;
          e[10] = z*z*nc +  c;
          e[11] = 0;

          e[12] = 0;
          e[13] = 0;
          e[14] = 0;
          e[15] = 1;
        }

        return this;
    }

    public Matrix4 SetTranslate(float x, float y, float z)
    {
        float[] e = Elements;

        e[0] = 1;  e[4] = 0;  e[8]  = 0;  e[12] = x;
        e[1] = 0;  e[5] = 1;  e[9]  = 0;  e[13] = y;
        e[2] = 0;  e[6] = 0;  e[10] = 1;  e[14] = z;
        e[3] = 0;  e[7] = 0;  e[11] = 0;  e[15] = 1;

        return this;
    }

    public Matrix4 Translate(float x, float y, float z)
    {
        float[] e = Elements;
        
        e[12] += e[0] * x + e[4] * y + e[8]  * z;
        e[13] += e[1] * x + e[5] * y + e[9]  * z;
        e[14] += e[2] * x + e[6] * y + e[10] * z;
        e[15] += e[3] * x + e[7] * y + e[11] * z;

        return this;
    }

    /**
     * Transpose the matrix.
     * @return this
     */
    public Matrix4 Transpose()
    {
        float t;
        
        float[] e = Elements;
        
        t = e[ 1];  e[ 1] = e[ 4];  e[ 4] = t;
        t = e[ 2];  e[ 2] = e[ 8];  e[ 8] = t;
        t = e[ 3];  e[ 3] = e[12];  e[12] = t;
        t = e[ 6];  e[ 6] = e[ 9];  e[ 9] = t;
        t = e[ 7];  e[ 7] = e[13];  e[13] = t;
        t = e[11];  e[11] = e[14];  e[14] = t;
        
        return this;
    }
}
