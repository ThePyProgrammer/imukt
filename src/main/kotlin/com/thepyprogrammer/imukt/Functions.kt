package com.thepyprogrammer.imukt

import kotlin.math.sqrt

/**
 * Get a 1D x-sized Zero Double Array
 */
fun zeros(x: Int) = run {
    val newList = mutableListOf<Double>()
    for (i in 0..x) {
        newList.add(0.0)
    }
    newList.toTypedArray()
}

/**
 * Get a 2D row x col Zero Double Array
 */
fun zeros(row: Int, col: Int): Array<Array<Double>> = run {
    val newList = mutableListOf<Array<Double>>()
    for (i in 0..row) {
        newList.add(zeros(col))
    }
    newList.toTypedArray()
}

fun invSqrt(x: Float): Float {
//    val halfx = 0.5f * x;
//    var y = x;
//    var i = y.toLong()
//    i = 0x5f3759df - (i shr 1);
//    y = i.toFloat()
//    y *= 1.5f - (halfx * y * y)
//    y *= 1.5f - (halfx * y * y)
//    return y
    return 1f / sqrt(x)
}