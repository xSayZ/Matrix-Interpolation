using System;
using System.Collections.Generic;
using UnityEngine;
using Vectors;

[ExecuteAlways, RequireComponent(typeof(VectorRenderer))]
public class MatrixInterpolation : MonoBehaviour
{
    private VectorRenderer vectors;

    private List<Vector3> verts;
    private List<Vector3> originalVerts;

    [SerializeField, HideInInspector] public Matrix4x4 startMatrix;
    [SerializeField, HideInInspector] public Matrix4x4 endMatrix;
    [SerializeField, HideInInspector] public Matrix4x4 changeMatrix;

    [SerializeField][Range(0, 1.0f)] private float interpolationFactor;

    [Header("Toggle Interpolation")]
    [SerializeField] private bool interpolateRotation;
    [SerializeField] private bool interpolateScale;
    [SerializeField] private bool interpolateTranslation;

    // Start is called before the first frame update
    void Start() {
        InitializeVectors();
        InitializeOriginalVerts();
    }

    void InitializeOriginalVerts() {
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

    void InitializeVectors() {
        if (!TryGetComponent<VectorRenderer>(out vectors))
        {
            vectors = gameObject.AddComponent<VectorRenderer>();
        }
    }

    // Update is called once per frame
    void Update()
    {
        Matrix4x4 rotationMatrix = CalculateRotationMatrix();
        Matrix4x4 scaleMatrix = CalculateScaleMatrix(rotationMatrix);

        // Transpose the changeMatrix
        changeMatrix = scaleMatrix;
        changeMatrix = changeMatrix.Transpose();
        ClearTranslationComponents();

        // Check if interpolation of translation is enabled
        if (interpolateTranslation) {
            InterpolateTranslation();
        }
        else {
            UseStartMatrixTranslation();
        }

        // Multiply the original vertices of the cube with the changeMatrix to get the transformed vertices
        verts = MultiplyVertsMatrix(changeMatrix, originalVerts);

        // Draw the transformed cube
        DrawCube();
        
    }

    void UseStartMatrixTranslation() {
        // If interpolation of translation is disabled, use the startMatrix's translation directly
        changeMatrix.m03 = startMatrix.m03;
        changeMatrix.m13 = startMatrix.m13;
        changeMatrix.m23 = startMatrix.m23;
    }

    void InterpolateTranslation() {
        // Interpolate the translation values and update the changeMatrix
        Vector3 lerpResult = MatrixOperations.Lerp( new Vector3(startMatrix.m03, startMatrix.m13, startMatrix.m23),  new Vector3(endMatrix.m03, endMatrix.m13, endMatrix.m23),interpolationFactor);
        changeMatrix.m03 = lerpResult.x;
        changeMatrix.m13 = lerpResult.y;
        changeMatrix.m23 = lerpResult.z;
    }

    void ClearTranslationComponents() {
        // Clear the translation components of the changeMatrix
        Vector4 mRow0 = changeMatrix.GetColumn(0);
        Vector4 mRow1 = changeMatrix.GetColumn(1);
        Vector4 mRow2 = changeMatrix.GetColumn(2);
        mRow0.w = 0; mRow1.w = 0; mRow2.w = 0;
        changeMatrix.SetColumn(0, mRow0);
        changeMatrix.SetColumn(1, mRow1);
        changeMatrix.SetColumn(2, mRow2);
    }

    Matrix4x4 CalculateRotationMatrix() {
        changeMatrix = startMatrix;

        // Calculate the rotationMatrix using Quaternion interpolation
        Matrix4x4 rotationMatrix = MatrixOperations.QuaternionToMatrix(MatrixOperations.MatrixToQuaternion(startMatrix));

        // Check if interpolation of rotation is enabled
        if (interpolateRotation)
        {
            // Calculate Quaternion values for startMatrix and endMatrix
            Quaternion q = MatrixOperations.MatrixToQuaternion(startMatrix);
            Quaternion b = MatrixOperations.MatrixToQuaternion(endMatrix);

            // Interpolate the rotation using Slerp and update the rotationMatrix
            rotationMatrix = MatrixOperations.QuaternionToMatrix(MatrixOperations.Slerp(q, b, interpolationFactor));
        }
        else
        {
            // If interpolation of rotation is disabled, use the startMatrix's rotation directly
            Quaternion q = MatrixOperations.MatrixToQuaternion(startMatrix);
            rotationMatrix = MatrixOperations.QuaternionToMatrix(q);
        }

        // Invert the rotationMatrix
        rotationMatrix = rotationMatrix.inverse;
        return rotationMatrix;
    }
    
    Matrix4x4 CalculateScaleMatrix(Matrix4x4 rotationMatrix) {
        // Calculate the scaleMatrix for the startMatrix
        float scaleXstart = startMatrix.GetColumn(0).magnitude;
        float scaleYstart = startMatrix.GetColumn(1).magnitude;
        float scaleZstart = startMatrix.GetColumn(2).magnitude;
        Matrix4x4 scaleMatrix = MatrixOperations.ScaleMatrix4X4(scaleXstart, scaleYstart,scaleZstart);

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
            float scaleX = MatrixOperations.Lerp(scaleXstart, scaleXend, interpolationFactor);
            float scaleY = MatrixOperations.Lerp(scaleYstart, scaleYend, interpolationFactor);
            float scaleZ = MatrixOperations.Lerp(scaleZstart, scaleZend, interpolationFactor);
            scaleMatrix = MatrixOperations.ScaleMatrix4X4(scaleX, scaleY, scaleZ);
        }

        // Combine the scaleMatrix and rotationMatrix
        scaleMatrix *= rotationMatrix;
        return scaleMatrix;
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