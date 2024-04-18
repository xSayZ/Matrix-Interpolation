using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class Matrix 
{
    public static Matrix4x4 Transpose(this Matrix4x4 m)
    {
        Vector4 c0 = m.GetColumn(0);
        Vector4 c1 = m.GetColumn(1);
        Vector4 c2 = m.GetColumn(2);
        
        m.SetRow(0, c0);
        m.SetRow(1, c1);
        m.SetRow(2, c2);
        return m;
    }
}