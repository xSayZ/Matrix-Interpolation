using UnityEditor;
using UnityEngine;

[CustomEditor((typeof(MatrixInterpolation)))]
public class MatrixInterpolationEditor : Editor
{
    void OnSceneGUI()
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
    void DrawMatrixProperty(string label, ref Matrix4x4 matrix)
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