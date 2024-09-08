using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class MatrixOperations
{
    public static Matrix4x4 ScaleMatrix4X4(float scaleX, float scaleY, float scaleZ)
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
    public static Quaternion Slerp(Quaternion q0, Quaternion q1, float t)
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
    public static Quaternion ScalarQuaternionMultiplication(Quaternion q, float scalar)
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
    public static Vector3 Lerp(Vector3 a, Vector3 b, float t)
    {
        // Perform linear interpolation between vectors 'a' and 'b' based on 't'
        // The resulting vector is (1 - t) times vector 'a' plus 't' times vector 'b'
        return (1 - t) * a + t * b;
    }

    // Linear interpolation method for float values
    public static float Lerp(float a, float b, float t)
    {
        // Perform linear interpolation between 'a' and 'b' based on 't'
        // The result is (1 - t) times 'a' plus 't' times 'b'
        return (1 - t) * a + t * b;
    }
    
    public static Quaternion MatrixToQuaternion(Matrix4x4 m)
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

    public static Matrix4x4 QuaternionToMatrix(Quaternion q)
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
    

    public static Vector3 MatrixToScale(Matrix4x4 m)
    {
        Vector3 scaleV = new Vector3();
        scaleV.x = m.GetColumn(0).magnitude;
        scaleV.y = m.GetColumn(1).magnitude;
        scaleV.z = m.GetColumn(2).magnitude;
        return scaleV;
    }
}
