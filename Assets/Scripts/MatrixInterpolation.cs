using System;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;
using Vectors;

[ExecuteAlways, RequireComponent(typeof(VectorRenderer))]
public class MatrixInterpolation : MonoBehaviour
{
    private VectorRenderer vectors;

    private List<Vector3> verts;
    private List<Vector3> originalVerts;

    [SerializeField, HideInInspector] internal Vector3 startMatrixTarget;
    [SerializeField, HideInInspector] internal Matrix4x4 startMatrix;
    [SerializeField, HideInInspector] internal Matrix4x4 endMatrix;
    [SerializeField, HideInInspector] internal Matrix4x4 changeMatrix;

    [SerializeField][Range(0, 1.0f)] private float interpolationFactor;

    [Header("Toggle Interpolation")]
    [SerializeField] private bool interpolateRotation;
    [SerializeField] private bool interpolateScale;
    [SerializeField] private bool interpolateTranslation;

    // Start is called before the first frame update
    void Start()
    {
        if (!TryGetComponent<VectorRenderer>(out vectors))
        {
            vectors = gameObject.AddComponent<VectorRenderer>();
        }

        originalVerts = new List<Vector3>();
        
        // Defining the 1st cube vertices
        originalVerts.Add(new Vector3(-0.5f, -0.5f, -0.5f));
        originalVerts.Add(new Vector3(-0.5f, -0.5f, 0.5f));
        originalVerts.Add(new Vector3(0.5f, -0.5f, -0.5f));
        originalVerts.Add(new Vector3(0.5f, -0.5f, 0.5f));
        originalVerts.Add(new Vector3(-0.5f, 0.5f, -0.5f));
        originalVerts.Add(new Vector3(-0.5f, 0.5f, 0.5f));
        originalVerts.Add(new Vector3(0.5f, 0.5f, -0.5f));
        originalVerts.Add(new Vector3(0.5f, 0.5f, 0.5f));
        originalVerts.Add(new Vector3(0,0,0));
        originalVerts.Add(new Vector3(1,0,0));
        originalVerts.Add(new Vector3(0,1,0));
        originalVerts.Add(new Vector3(0,0,1));
    }

    // Update is called once per frame
    void Update()
    {
        changeMatrix = startMatrix;

        // Calculate the rotationMatrix using Quaternion interpolation
        Matrix4x4 rotationMatrix = MatrixInterpolationEditor.QuaternionToMatrix(MatrixInterpolationEditor.MatrixToQuaternion(startMatrix));

        // Check if interpolation of rotation is enabled
        if (interpolateRotation)
        {
            // Calculate Quaternion values for startMatrix and endMatrix
            Quaternion q = MatrixInterpolationEditor.MatrixToQuaternion(startMatrix);
            Quaternion b = MatrixInterpolationEditor.MatrixToQuaternion(endMatrix);

            // Interpolate the rotation using Slerp and update the rotationMatrix
            rotationMatrix = MatrixInterpolationEditor.QuaternionToMatrix(Slerp(q, b, interpolationFactor));
        }
        else
        {
            // If interpolation of rotation is disabled, use the startMatrix's rotation directly
            Quaternion q = MatrixInterpolationEditor.MatrixToQuaternion(startMatrix);
            rotationMatrix = MatrixInterpolationEditor.QuaternionToMatrix(q);
        }

        // Invert the rotationMatrix
        rotationMatrix = rotationMatrix.inverse;

        // Calculate the scaleMatrix for the startMatrix
        float scaleXstart = startMatrix.GetColumn(0).magnitude;
        float scaleYstart = startMatrix.GetColumn(1).magnitude;
        float scaleZstart = startMatrix.GetColumn(2).magnitude;
        Matrix4x4 scaleMatrix = ScaleMatrix4X4(scaleXstart, scaleYstart,scaleZstart);

        // Check if interpolation of scale is enabled
        if (interpolateScale)
        {
            // Calculate scale values for both startMatrix and endMatrix
            scaleXstart = startMatrix.GetColumn(0).magnitude;
            scaleYstart = startMatrix.GetColumn(1).magnitude;
            scaleZstart = startMatrix.GetColumn(2).magnitude;
            float scaleXend = endMatrix.GetColumn(0).magnitude;
            float scaleYend = endMatrix.GetColumn(1).magnitude;
            float scaleZend = endMatrix.GetColumn(2).magnitude;

            // Interpolate the scale using Lerp and update the scaleMatrix
            float scaleX = Lerp(scaleXstart, scaleXend, interpolationFactor);
            float scaleY = Lerp(scaleYstart, scaleYend, interpolationFactor);
            float scaleZ = Lerp(scaleZstart, scaleZend, interpolationFactor);
            scaleMatrix = ScaleMatrix4X4(scaleX, scaleY, scaleZ);
        }

        // Combine the scaleMatrix and rotationMatrix
        scaleMatrix *= rotationMatrix;

        // Transpose the changeMatrix
        changeMatrix = scaleMatrix;
        changeMatrix = changeMatrix.Transpose();

        // Clear the translation components of the changeMatrix
        Vector4 mRow0 = changeMatrix.GetColumn(0);
        Vector4 mRow1 = changeMatrix.GetColumn(1);
        Vector4 mRow2 = changeMatrix.GetColumn(2);
        mRow0.w = 0; mRow1.w = 0; mRow2.w = 0;
        changeMatrix.SetColumn(0, mRow0);
        changeMatrix.SetColumn(1, mRow1);
        changeMatrix.SetColumn(2, mRow2);

        // Check if interpolation of translation is enabled
        if (interpolateTranslation)
        {
            // Interpolate the translation values and update the changeMatrix
            Vector3 lerpResult = Lerp( new Vector3(startMatrix.m03, startMatrix.m13, startMatrix.m23),  new Vector3(endMatrix.m03, endMatrix.m13, endMatrix.m23),interpolationFactor);
            changeMatrix.m03 = lerpResult.x;
            changeMatrix.m13 = lerpResult.y;
            changeMatrix.m23 = lerpResult.z;
        }
        else
        {
            // If interpolation of translation is disabled, use the startMatrix's translation directly
            changeMatrix.m03 = startMatrix.m03;
            changeMatrix.m13 = startMatrix.m13;
            changeMatrix.m23 = startMatrix.m23;
        }

        // Multiply the original vertices of the cube with the changeMatrix to get the transformed vertices
        verts = MultiplyVertsMatrix(changeMatrix, originalVerts);

        // Draw the transformed cube
        DrawCube();
        
    }

    // Method to create a scaling Matrix4x4
    private Matrix4x4 ScaleMatrix4X4(float scaleX, float scaleY, float scaleZ)
    {
        // Create a new identity Matrix4x4, representing no transformation initially
        Matrix4x4 m = Matrix4x4.identity;

        // Set the scaling factors on the diagonal elements of the Matrix4x4
        // The diagonal elements represent scaling along the x, y, and z axes
        m.m00 = scaleX; // Scaling factor for the x-axis
        m.m11 = scaleY; // Scaling factor for the y-axis
        m.m22 = scaleZ; // Scaling factor for the z-axis

        // Return the resulting scaling Matrix4x4
        return m;
    }

    // Method to perform spherical linear interpolation (SLERP) between two Quaternions
    private Quaternion Slerp(Quaternion q0, Quaternion q1, float t)
    {
        // Create a new Quaternion to store the result of the interpolation
        Quaternion finalQ = new Quaternion();

        // Calculate the dot product between the two Quaternions
        float cosW = Quaternion.Dot(q0, q1);

        // Check if the angle between the two Quaternions is greater than 90 degrees (cosine < 0)
        if (cosW < 0.0f)
        {
            // If the angle is greater than 90 degrees, negate one of the Quaternions
            // to ensure the shortest path for the interpolation
            ScalarQuaternionMultiplication(q1, -1);
        }

        // Variables to store the interpolation weights for q0 and q1
        float k0, k1;

        // Check if the angle between the two Quaternions is very close to 0 degrees (cosine ~ 1)
        if (cosW > 0.9999f)
        {
            // If the angle is very close to 0 degrees, perform a simple linear interpolation
            k0 = 1.0f - t;
            k1 = t;
        }
        else
        {
            // Calculate the angle (omega) and sine of the angle (sinOmega) between the Quaternions
            float sinOmega = Mathf.Sqrt(1.0f - cosW *cosW);

            // Calculate the reciprocal of sinOmega
            float omega = Mathf.Atan2(sinOmega, cosW);
            float inverseSinOmega = 1.0f / sinOmega;

            // Calculate the interpolation weights using spherical linear interpolation formula
            k0 = Mathf.Sin((1.0f - t) * omega) * inverseSinOmega;
            k1 = Mathf.Sin(t * omega) * inverseSinOmega;
        }

        // Interpolate the components of the Quaternions and store the result in finalQ
        finalQ.x = q0.x * k0 + q1.x * k1;
        finalQ.y = q0.y * k0 + q1.y * k1;
        finalQ.z = q0.z * k0 + q1.z * k1;
        finalQ.w = q0.w * k0 + q1.w * k1;

        // Return the resulting interpolated Quaternion
        return finalQ;
    }

    // Method to perform scalar multiplication on a Quaternion
    private Quaternion ScalarQuaternionMultiplication(Quaternion q, float scalar)
    {
        // Scale the individual components of the Quaternion by the scalar value
        q.x *= scalar;
        q.y *= scalar;
        q.z *= scalar;
        q.w *= scalar;

        // Return the resulting scaled Quaternion
        return q;
    }

    // Linear interpolation method for Vector3
    private Vector3 Lerp(Vector3 a, Vector3 b, float t)
    {
        // Perform linear interpolation between vectors 'a' and 'b' based on 't'
        // The resulting vector is (1 - t) times vector 'a' plus 't' times vector 'b'
        return (1 - t) * a + t * b;
    }

    // Linear interpolation method for float values
    private float Lerp(float a, float b, float t)
    {
        // Perform linear interpolation between 'a' and 'b' based on 't'
        // The result is (1 - t) times 'a' plus 't' times 'b'
        return (1 - t) * a + t * b;
    }
    private void DrawCube()
    {
        using (vectors.Begin())
        {
            vectors.Draw(verts[0], verts[1], Color.red);
            vectors.Draw(verts[2], verts[3], Color.red);
            vectors.Draw(verts[0], verts[2], Color.yellow);
            vectors.Draw(verts[1], verts[3], Color.yellow);
            vectors.Draw(verts[0], verts[4], Color.blue);
            vectors.Draw(verts[1], verts[5], Color.blue);
            vectors.Draw(verts[2], verts[6], Color.blue);
            vectors.Draw(verts[3], verts[7], Color.blue);
            vectors.Draw(verts[4], verts[5], Color.red);
            vectors.Draw(verts[6], verts[7], Color.red);
            vectors.Draw(verts[4], verts[6], Color.yellow);
            vectors.Draw(verts[5], verts[7], Color.yellow);
            vectors.Draw(verts[8], verts[9], Color.red);
            vectors.Draw(verts[8],verts[10], Color.green);
            vectors.Draw(verts[8], verts[11], Color.blue);
        }
    }

    // Method to transform a list of Vector3 vertices using a given Matrix4x4
    private List<Vector3> MultiplyVertsMatrix(Matrix4x4 mat, List<Vector3> vertList)
    {
        // Create a new list to store the transformed vertices
        List<Vector3> result = new List<Vector3>();

        // Iterate through each vertex in the original vertex list
        for (int i = 0; i < vertList.Count; i++)
        {
            // Use the Matrix4x4 method MultiplyPoint to transform each vertex
            // The transformed vertex is added to the 'result' list
            result.Add(mat.MultiplyPoint(vertList[i]));
        }

        // Return the list of transformed vertices
        return result;
    }
}

[CustomEditor((typeof(MatrixInterpolation)))]
public class MatrixInterpolationEditor : Editor
{
    static public Quaternion MatrixToQuaternion(Matrix4x4 m)
    {
       
        Vector4 mRow0 = m.GetColumn(0);
        Vector4 mRow1 = m.GetColumn(1);
        Vector4 mRow2 = m.GetColumn(2);
        
        m.SetColumn(0, mRow0.normalized);
        m.SetColumn(1, mRow1.normalized);
        m.SetColumn(2, mRow2.normalized);
        
        // Determine which of w, x, y, or z has the largest absolute value
        float traceW = m.m00 + m.m11 + m.m22;
        float traceX = m.m00 - m.m11 - m.m22;
        float traceY = m.m11 - m.m00 - m.m22;
        float traceZ = m.m22 - m.m00 - m.m11;
        float w = 0, x = 0, y = 0, z = 0;
        int biggestIndex = 0;
        float traceBiggest = traceW;
        if (traceX > traceBiggest) {
            traceBiggest = traceX;
            biggestIndex = 1;
        }
        if (traceY > traceBiggest) {
            traceBiggest = traceY;
            biggestIndex = 2;
        }
        if (traceZ > traceBiggest) {
            traceBiggest = traceZ;
            biggestIndex = 3;
        }

        // Perform square root and division
        float sqrtTraceBiggest = Mathf.Sqrt(traceBiggest + 1.0f) * 0.5f;
        float fraction = 0.25f / sqrtTraceBiggest;
        
        // Apply table to compute quaternion values
        switch (biggestIndex) {
            case 0:
                w = sqrtTraceBiggest;
                x = (m.m21 - m.m12) * fraction;
                y = (m.m02 - m.m20) * fraction;
                z = (m.m10 - m.m01) * fraction;
                break;

            case 1:
                x = sqrtTraceBiggest;
                w = (m.m21 - m.m12) * fraction;
                y = (m.m10 + m.m01) * fraction;
                z = (m.m02 + m.m20) * fraction;
                break;

            case 2:
                y = sqrtTraceBiggest;
                w = (m.m02 - m.m20) * fraction;
                x = (m.m10 + m.m01) * fraction;
                z = (m.m21 + m.m12) * fraction;
                break;

            case 3:
                z = sqrtTraceBiggest;
                w = (m.m10 - m.m01) * fraction;
                x = (m.m02 + m.m20) * fraction;
                y = (m.m21 + m.m12) * fraction; 
                break;
        }

        return new Quaternion(x, y, z, w);
    }

    static public Matrix4x4 QuaternionToMatrix(Quaternion q)
    {
        Matrix4x4 m = Matrix4x4.identity;

        float x = q.x; float y = q.y;
        float z = q.z; float w = q.w;

        float x2 = x * 2;
        float y2 = y * 2;
        float z2 = z * 2;
        float w2 = w * 2;

        m.m00 = 1 - y2 * y - z2 * z;
        m.m10 = x2 * y - w2 * z;
        m.m20 = x2 * z + w2 * y;

        m.m01 = x2 * y + w2 * z;
        m.m11 = 1 - x2 * x - z2 * z;
        m.m21 = y2 * z - w2 * x;

        m.m02 = x2 * z - w2 * y;
        m.m12 = y2 * z + w2 * x;
        m.m22 = 1 - x2 * x - y2 * y;

        return m.Transpose();
    }

    private Matrix4x4 ScaleMatrix4X4(Vector3 v)
    {
        Matrix4x4 m = Matrix4x4.identity;
        m.m00 = v.x;
        m.m11 = v.y;
        m.m22 = v.z;
        return m;
    }

    static public Vector3 MatrixToScale(Matrix4x4 m)
    {
        Vector3 scaleV = new Vector3();
        scaleV.x = m.GetColumn(0).magnitude;
        scaleV.y = m.GetColumn(1).magnitude;
        scaleV.z = m.GetColumn(2).magnitude;
        return scaleV;
    }
    private void OnSceneGUI()
    {
        var script = target as MatrixInterpolation;
        if (!script) return;

        EditorGUI.BeginChangeCheck();
        var newTargetStart = Handles.PositionHandle(script.startMatrix.GetTranslation(), script.transform.rotation);
        var newScaleStart = Handles.ScaleHandle(script.startMatrix.GetScale(), script.startMatrix.GetTranslation(), script.startMatrix.GetRotation(), 0.5f);
        var newRotationStart = Handles.RotationHandle(script.startMatrix.GetRotation(), script.startMatrix.GetTranslation());

        if (EditorGUI.EndChangeCheck())
        {
            Undo.RecordObject(target, "Moved Start Matrix");
            script.startMatrix = Matrix4x4.TRS(newTargetStart, newRotationStart, newScaleStart);
            EditorUtility.SetDirty(script);
        }

        EditorGUI.BeginChangeCheck();
        var newTargetEnd = Handles.PositionHandle(script.endMatrix.GetTranslation(), script.transform.rotation);
        var newScaleEnd = Handles.ScaleHandle(script.endMatrix.GetScale(), script.endMatrix.GetTranslation(), script.endMatrix.GetRotation(), 0.5f);
        var newRotationEnd = Handles.RotationHandle(script.endMatrix.GetRotation(), script.endMatrix.GetTranslation());

        if (EditorGUI.EndChangeCheck())
        {
            Undo.RecordObject(target, "Moved End Matrix");
            script.endMatrix = Matrix4x4.TRS(newTargetEnd, newRotationEnd, newScaleEnd);
            EditorUtility.SetDirty(script);
        }
    }

    public override void OnInspectorGUI()
    {
        // Call the base class implementation to draw the default inspector GUI.
        base.OnInspectorGUI();

        // Get the target component as MatrixInterpolation.
        var script = target as MatrixInterpolation;

        // If the target is not a MatrixInterpolation component, return early and do nothing.
        if (!script) return;

        // Draw the editable properties for the start and end matrices.
        DrawMatrixProperty("Start Matrix", ref script.startMatrix);
        DrawMatrixProperty("End Matrix", ref script.endMatrix);

        // Add a space to visually separate matrix properties from determinant properties.
        EditorGUILayout.Space();

        // Display the determinants of the start, change, and end matrices.
        EditorGUILayout.PrefixLabel("Start Matrix Determinant");
        EditorGUILayout.FloatField(script.startMatrix.determinant);
        EditorGUILayout.PrefixLabel("Change Matrix Determinant");
        EditorGUILayout.FloatField(script.changeMatrix.determinant);
        EditorGUILayout.PrefixLabel("End Matrix Determinant");
        EditorGUILayout.FloatField(script.endMatrix.determinant);
    }

    // Helper method to draw the editable matrix property in the inspector.
    private void DrawMatrixProperty(string label, ref Matrix4x4 matrix)
    {
        // Begin checking for changes to the matrix property.
        EditorGUI.BeginChangeCheck();

        // Draw the label and the matrix fields in a horizontal layout.
        EditorGUILayout.BeginHorizontal();
        EditorGUILayout.PrefixLabel(label);
        EditorGUILayout.BeginVertical();

        // Create a temporary matrix to hold the modified values.
        var result = matrix;

        // Loop to draw the individual elements of the matrix.
        for (int i = 0; i < 4; i++)
        {
            EditorGUILayout.BeginHorizontal();
            for (int j = 0; j < 4; j++)
            {
                // Draw the float field for each matrix element and update the temporary matrix.
                result[i, j] = EditorGUILayout.FloatField(matrix[i, j]);
            }
            EditorGUILayout.EndHorizontal();
        }

        EditorGUILayout.EndVertical();
        EditorGUILayout.EndHorizontal();

        // If any changes were detected during the loop, record the object for Undo support.
        if (EditorGUI.EndChangeCheck())
        {
            Undo.RecordObject(target, "Change Matrix");

            // Update the original matrix with the modified values.
            matrix = result;

            // Mark the target object as dirty so changes are saved.
            EditorUtility.SetDirty(target);
        }
    }

}

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