using UnityEngine;

public static class Matrix4x4Extensions
{
    public static Vector3 GetTranslation(this Matrix4x4 matrix)
    {
        return matrix.GetColumn(3);
    }

    public static Vector3 GetScale(this Matrix4x4 matrix)
    {
        Vector3 scale;
        scale.x = new Vector4(matrix.m00, matrix.m10, matrix.m20, matrix.m30).magnitude;
        scale.y = new Vector4(matrix.m01, matrix.m11, matrix.m21, matrix.m31).magnitude;
        scale.z = new Vector4(matrix.m02, matrix.m12, matrix.m22, matrix.m32).magnitude;
        return scale;
    }

    public static Quaternion GetRotation(this Matrix4x4 matrix)
    {
        // Extract rotation from the matrix using LookRotation method
        // Normalize the result as LookRotation can return a non-normalized quaternion
        return Quaternion.LookRotation(matrix.GetColumn(2), matrix.GetColumn(1)).normalized;
    }
}