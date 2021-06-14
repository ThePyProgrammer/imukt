package com.thepyprogrammer.imukt

import kotlin.math.asin
import kotlin.math.atan2
import kotlin.math.sqrt

// Source: https://github.com/arduino-libraries/MadgwickAHRS
class Madgwick {
    var beta = 0.1f

    var q0 = 1f
    var q1 = 0f
    var q2 = 0f
    var q3 = 0f

    var invSampleFreq = 1f / 512f

    private var roll = 0f
    private var pitch = 0f
    private var yaw: Float = 0f
    var anglesComputed = 0

    init {
        computeAngles()
    }

    fun begin(sampleFrequency: Float) {
        invSampleFreq = 1.0f / sampleFrequency
    }

    fun update(gxInDeg: Float, gyInDeg: Float, gzInDeg: Float, axVal: Float, ayVal: Float, azVal: Float, mxVal: Float, myVal: Float, mzVal: Float) {
        if(mxVal == 0f && myVal == 0f && mzVal == 0f) {
            updateIMU(gxInDeg, gyInDeg, gzInDeg, axVal, ayVal, azVal)
        }

        // Convert gyroscope degrees/sec to radians/sec
        val gx = gxInDeg * 0.0174533f
        val gy = gyInDeg * 0.0174533f
        val gz = gzInDeg * 0.0174533f

        var ax = axVal
        var ay = ayVal
        var az = azVal

        var mx = mxVal
        var my = myVal
        var mz = mzVal

        // Rate of change of quaternion from gyroscope
        var qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz)
        var qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy)
        var qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx)
        var qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx)

        var recipNorm: Float

        // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
        if(!((ax == 0f) && (ay == 0f) && (az == 0f))) {

            // Normalise accelerometer measurement
            recipNorm = invSqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm

            // Normalise magnetometer measurement
            recipNorm = invSqrt(mx * mx + my * my + mz * mz)
            mx *= recipNorm
            my *= recipNorm
            mz *= recipNorm

            // Auxiliary variables to avoid repeated arithmetic
            val _2q0mx = 2.0f * q0 * mx
            val _2q0my = 2.0f * q0 * my
            val _2q0mz = 2.0f * q0 * mz
            val _2q1mx = 2.0f * q1 * mx
            val _2q0 = 2.0f * q0
            val _2q1 = 2.0f * q1
            val _2q2 = 2.0f * q2
            val _2q3 = 2.0f * q3
            val _2q0q2 = 2.0f * q0 * q2
            val _2q2q3 = 2.0f * q2 * q3
            val q0q0 = q0 * q0
            val q0q1 = q0 * q1
            val q0q2 = q0 * q2
            val q0q3 = q0 * q3
            val q1q1 = q1 * q1
            val q1q2 = q1 * q2
            val q1q3 = q1 * q3
            val q2q2 = q2 * q2
            val q2q3 = q2 * q3
            val q3q3 = q3 * q3

            // Reference direction of Earth's magnetic field
            val hx = mx * q0q0 - _2q0my * q3 + _2q0mz * q2 + mx * q1q1 + _2q1 * my * q2 + _2q1 * mz * q3 - mx * q2q2 - mx * q3q3
            val hy = _2q0mx * q3 + my * q0q0 - _2q0mz * q1 + _2q1mx * q2 - my * q1q1 + my * q2q2 + _2q2 * mz * q3 - my * q3q3
            val _2bx = sqrt(hx * hx + hy * hy)
            val _2bz = -_2q0mx * q2 + _2q0my * q1 + mz * q0q0 + _2q1mx * q3 - mz * q1q1 + _2q2 * my * q3 - mz * q2q2 + mz * q3q3
            val _4bx = 2.0f * _2bx
            val _4bz = 2.0f * _2bz

            // Gradient decent algorithm corrective step
            var s0: Float = -_2q2 * (2.0f * q1q3 - _2q0q2 - ax) + _2q1 * (2.0f * q0q1 + _2q2q3 - ay) - _2bz * q2 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q3 + _2bz * q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz)
            var s1: Float = _2q3 * (2.0f * q1q3 - _2q0q2 - ax) + _2q0 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q1 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + _2bz * q3 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q2 + _2bz * q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q3 - _4bz * q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz)
            var s2: Float = -_2q0 * (2.0f * q1q3 - _2q0q2 - ax) + _2q3 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q2 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + (-_4bx * q2 - _2bz * q0) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q1 + _2bz * q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q0 - _4bz * q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz)
            var s3: Float = _2q1 * (2.0f * q1q3 - _2q0q2 - ax) + _2q2 * (2.0f * q0q1 + _2q2q3 - ay) + (-_4bx * q3 + _2bz * q1) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q0 + _2bz * q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz)
            recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3) // normalise step magnitude
            s0 *= recipNorm
            s1 *= recipNorm
            s2 *= recipNorm
            s3 *= recipNorm

            // Apply feedback step
            qDot1 -= beta * s0
            qDot2 -= beta * s1
            qDot3 -= beta * s2
            qDot4 -= beta * s3
        }

        // Integrate rate of change of quaternion to yield quaternion
        q0 += qDot1 * invSampleFreq
        q1 += qDot2 * invSampleFreq
        q2 += qDot3 * invSampleFreq
        q3 += qDot4 * invSampleFreq

        // Normalise quaternion
        recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
        q0 *= recipNorm
        q1 *= recipNorm
        q2 *= recipNorm
        q3 *= recipNorm
        anglesComputed = 0
    }

    fun updateIMU(gxInDeg: Float, gyInDeg: Float, gzInDeg: Float, axVal: Float, ayVal: Float, azVal: Float) {
        // Convert gyroscope degrees/sec to radians/sec
        val gx = gxInDeg * 0.0174533f
        val gy = gyInDeg * 0.0174533f
        val gz = gzInDeg * 0.0174533f

        var ax = axVal
        var ay = ayVal
        var az = azVal

        var recipNorm: Float

        // Rate of change of quaternion from gyroscope
        var qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz)
        var qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy)
        var qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx)
        var qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx)

        // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
        if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

            // Normalise accelerometer measurement
            recipNorm = invSqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm

            // Auxiliary variables to avoid repeated arithmetic
            val _2q0 = 2.0f * q0
            val _2q1 = 2.0f * q1
            val _2q2 = 2.0f * q2
            val _2q3 = 2.0f * q3
            val _4q0 = 4.0f * q0
            val _4q1 = 4.0f * q1
            val _4q2 = 4.0f * q2
            val _8q1 = 8.0f * q1
            val _8q2 = 8.0f * q2
            val q0q0 = q0 * q0
            val q1q1 = q1 * q1
            val q2q2 = q2 * q2
            val q3q3 = q3 * q3

            // Gradient decent algorithm corrective step
            var s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay
            var s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az
            var s2 = 4.0f * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az
            var s3 = 4.0f * q1q1 * q3 - _2q1 * ax + 4.0f * q2q2 * q3 - _2q2 * ay
            recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3) // normalise step magnitude
            s0 *= recipNorm
            s1 *= recipNorm
            s2 *= recipNorm
            s3 *= recipNorm

            // Apply feedback step
            qDot1 -= beta * s0
            qDot2 -= beta * s1
            qDot3 -= beta * s2
            qDot4 -= beta * s3
        }

        // Integrate rate of change of quaternion to yield quaternion
        q0 += qDot1 * invSampleFreq
        q1 += qDot2 * invSampleFreq
        q2 += qDot3 * invSampleFreq
        q3 += qDot4 * invSampleFreq

        // Normalise quaternion
        recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
        q0 *= recipNorm
        q1 *= recipNorm
        q2 *= recipNorm
        q3 *= recipNorm
        anglesComputed = 0

    }

    fun computeAngles() {
        roll = atan2(q0*q1 + q2*q3, 0.5f - q1*q1 - q2*q2)
        pitch = asin(-2.0f * (q1*q3 - q0*q2))
        yaw = atan2(q1*q2 + q0*q3, 0.5f - q2*q2 - q3*q3)
        anglesComputed = 1
    }

    fun getRoll(): Float {
        if (anglesComputed == 0) computeAngles()
        return roll * 57.29578f
    }

    fun getPitch(): Float {
        if (anglesComputed == 0) computeAngles()
        return pitch * 57.29578f
    }

    fun getYaw(): Float {
        if (anglesComputed == 0) computeAngles()
        return yaw * 57.29578f + 180.0f
    }

    fun getRollRadians(): Float {
        if (anglesComputed == 0) computeAngles()
        return roll
    }

    fun getPitchRadians(): Float {
        if (anglesComputed == 0) computeAngles()
        return pitch
    }

    fun getYawRadians(): Float {
        if (anglesComputed == 0) computeAngles()
        return yaw
    }

}

