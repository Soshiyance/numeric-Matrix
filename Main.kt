package processor

import java.math.RoundingMode
import java.text.DecimalFormat
import kotlin.math.pow
import kotlin.math.round

fun main() {

    while (true) {

        printMenu()
        print("Your choice: ")
        val order = readln()
        if (order == "0") break
        // Add matrices
        if (order == "1") {
            print("Enter size of first matrix: ")
            val (firstRaw, firstCol) = readln().split(" ").map { it.toInt() }
            println("Enter first matrix:")
            val firstMatrix = Matrix(firstRaw, firstCol)
            firstMatrix.matrix = firstMatrix.matrixCreator()
            print("Enter size of second matrix: ")
            val (secondRaw, secondCol) = readln().split(" ").map { it.toInt() }
            println("Enter second matrix:")
            val secondMatrix = Matrix(secondRaw, secondCol)
            secondMatrix.matrix = secondMatrix.matrixCreator()
            println("The result is:")
            if (firstMatrix.isAddable(secondMatrix)) {
                val result = firstMatrix.sumMatrix(secondMatrix)
                result.printMatrix()
                println()
            } else {
                println("The operation cannot be performed.")
                println()
            }
        }
        // Multiply matrix by a constant
        if (order == "2") {
            print("Enter size of first matrix: ")
            val (firstRaw, firstCol) = readln().split(" ").map { it.toInt() }
            println("Enter first matrix:")
            val firstMatrix = Matrix(firstRaw, firstCol)
            firstMatrix.matrix = firstMatrix.matrixCreator()
            print("Enter constant: ")
            val scale = readln().toDouble()
            val result = firstMatrix.matrix.multiNum(scale)
            println("The result is:")
            result.printMatrix()
        }
        // Multiply matrices
        if (order == "3") {
            print("Enter size of first matrix: ")
            val (firstRaw, firstCol) = readln().split(" ").map { it.toInt() }
            println("Enter first matrix:")
            val firstMatrix = Matrix(firstRaw, firstCol)
            firstMatrix.matrix = firstMatrix.matrixCreator()
            print("Enter size of second matrix: ")
            val (secondRaw, secondCol) = readln().split(" ").map { it.toInt() }
            println("Enter second matrix:")
            val secondMatrix = Matrix(secondRaw, secondCol)
            secondMatrix.matrix = secondMatrix.matrixCreator()
            if (firstMatrix.isMultiple(secondMatrix)) {
                val result = firstMatrix.multiMatrix(secondMatrix)
                println("The result is:")
                result.printMatrix()
                println()
            } else {
                println("The operation cannot be performed.")
                println()
            }

        }
        // Transpose matrix
        if (order == "4") {
            printTransposeMenu()
            print("Your choice: ")
            val order = readln().toInt()
            print("Enter matrix size: ")
            val (raw, col) = readln().split(" ").map { it.toInt() }
            println("Enter matrix:")
            val matrix = Matrix(raw, col)
            matrix.matrix = matrix.matrixCreator()
            val transpose = when (order) {
                1 -> matrix.transposeDiagonally()
                2 -> matrix.transposeAntiDiagonally()
                3 -> matrix.transposeVertically()
                4 -> matrix.transposeHorizontally()
                else -> return

            }
            println("The result is:")
            transpose.printMatrix()
            println()
        }
        if (order == "5") {
            print("Your choice: ")
            val (raw, col) = readln().split(" ").map { it.toInt() }
            println("Enter matrix:")
            val matrix = Matrix(raw, col)
            matrix.matrix = matrix.matrixCreator()
            val determinant = matrix.matrix.calculateDeterminant()
            println("The result is:")
            println(determinant)
        }
        if (order == "6") {
            print("Enter matrix size: ")
            val (raw, col) = readln().split(" ").map { it.toInt() }
            println("Enter matrix:")
            val matrix = Matrix(raw, col)
            matrix.matrix = matrix.matrixCreator()
            val determinant = matrix.matrix.calculateDeterminant()
            if (determinant != 0.0) {
                val adjMatrix = matrix.matrix.adjMatrix()
                val invertMatrix = adjMatrix.transpose().multiNum(1 / determinant)
                println("The result is:")
                invertMatrix.printMatrix()
                println()
            } else {
                println("This matrix doesn't have an inverse.")
                println()
            }

        }

    }

}

fun MutableList<MutableList<Double>>.multiNum(n: Double): MutableList<MutableList<Double>> {
    for (raw in 0 until this.size) {
        for (col in 0 until this[raw].size) {
            this[raw][col] = this[raw][col] * n
        }
    }
    return this
}


class Matrix(raw: Int, col: Int) {

    private var raw: Int = raw
    private var col: Int = col
    var matrix = zeroMatrix()


    private fun zeroMatrix(): MutableList<MutableList<Double>> {

        return MutableList(this.raw) { MutableList(this.col) { 0.0 } }
    }

    fun matrixCreator(): MutableList<MutableList<Double>> {
        return MutableList(raw) { readln().split(" ").map { it.toDouble() }.toMutableList() }
    }


    fun isAddable(matrix: Matrix): Boolean {
        return this.raw == matrix.raw && this.col == matrix.col
    }

    fun isMultiple(matrix: Matrix): Boolean {
        return this.col == matrix.raw
    }

    private fun isZeroFraction(element: Double): Boolean {
        return element - round(element) == 0.0
    }

    fun printMatrix() {

        matrix.forEach { raw ->
            raw.forEach { col ->
                if (isZeroFraction(col)) {
                    print("${col.toInt()} ")
                } else print("$col ")
            }
            println()
        }
    }

    private fun takeColumn(numOfCol: Int): MutableList<Double> {
        val column = MutableList<Double>(this.matrix.size) { 0.0 }
        for (raw in 0 until this.matrix.size) {
            for (col in 0 until this.matrix[raw].size) {
                if (col == numOfCol) column[raw] = this.matrix[raw][col]
            }
        }
        return column
    }

    private fun multiSingleRaw(matrix: Matrix, raw: Int, col: Int): Double {
        var sum = 0.0


        for (j in 0 until this.matrix[raw].size) {
            sum += this.matrix[raw][j] * matrix.takeColumn(col)[j]
        }

        return sum
    }

    fun multiMatrix(matrix: Matrix): MutableList<MutableList<Double>> {
        val multoResult = MutableList(this.matrix.size) { MutableList(matrix.matrix[0].size) { 0.0 } }
        for (raw in 0 until multoResult.size) {
            for (col in 0 until multoResult[raw].size) {
                multoResult[raw][col] = multiSingleRaw(matrix, raw, col)
            }
        }
        return multoResult
    }

    fun sumMatrix(matrix: Matrix): MutableList<MutableList<Double>> {

        val sum = zeroMatrix()
        for (raw in 0 until this.matrix.size) {
            for (col in 0 until this.matrix[raw].size) {
                sum[raw][col] = this.matrix[raw][col] + matrix.matrix[raw][col]
            }
        }
        return sum
    }

    fun transposeDiagonally(): MutableList<MutableList<Double>> {
        val transpose = MutableList(this.col) { MutableList(this.raw) { 0.0 } }
        for (raw in 0 until this.matrix.size) {
            for (col in 0 until this.matrix[raw].size) {
                transpose[col][raw] = this.matrix[raw][col]
            }
        }
        return transpose
    }

    fun transposeAntiDiagonally(): MutableList<MutableList<Double>> {
        val transpose = MutableList(this.col) { MutableList(this.raw) { 0.0 } }
        for (raw in 0 until this.matrix.size) {
            for (col in 0 until this.matrix[raw].size) {
                transpose[this.matrix[raw].size - col - 1][this.matrix.size - raw - 1] = this.matrix[raw][col]
            }
        }
        return transpose
    }

    fun transposeVertically(): MutableList<MutableList<Double>> {
        val transpose = MutableList(this.raw) { MutableList(this.col) { 0.0 } }
        for (raw in 0 until this.matrix.size) {
            for (col in 0 until this.matrix[raw].size) {
                transpose[raw][col] = this.matrix[raw][this.matrix[raw].size - col - 1]
            }
        }

        return transpose
    }

    fun transposeHorizontally(): MutableList<MutableList<Double>> {
        val transpose = MutableList(this.raw) { MutableList(this.col) { 0.0 } }
        for (i in 0 until this.matrix.size) {
            transpose[i] = this.matrix[this.matrix.size - i - 1]
        }
        return transpose
    }


}


fun zeroMatrix(n: Int, m: Int): MutableList<MutableList<Int>> {

    return MutableList(n) { MutableList(m) { 0 } }
}


fun isZeroFraction(element: Double): Boolean {
    return element - round(element) == 0.0
}

fun MutableList<MutableList<Double>>.printMatrix() {
    this.joinToString("\n") { row ->
        row.joinToString(" ") {"%6.2f".format(it)} }.also { print(it) }
    println()
}



fun Double.round(decimals: Int): Double {
    val df = DecimalFormat("0.00")
    df.roundingMode = RoundingMode.DOWN

    return df.format(this).toDouble()
}

fun printMenu() {
    println("1. Add matrices")
    println("2. Multiply matrix by a constant")
    println("3. Multiply matrices")
    println("4. Transpose matrix")
    println("5. Calculate a determinant")
    println("6. Inverse matrix")
    println("0. Exit")
}

fun printTransposeMenu() {
    println("1. Main diagonal")
    println("2. Side diagonal")
    println("3. Vertical line")
    println("4. Horizontal line")
}

fun MutableList<MutableList<Double>>.findSubMatrix(rawIndex: Int, colIndex: Int): MutableList<MutableList<Double>> {
    val subMatrix = MutableList(this.size - 1) { MutableList(this.size - 1) { 0.0 } }
    var subRaw = 0
    var subCol = 0
    raw@ for (raw in 0 until this.size) {
        subCol = 0
        for (col in 0 until this[raw].size) {
            if (col == colIndex) {
                continue
            } else if (raw == rawIndex) {
                continue@raw
            } else {
                subMatrix[subRaw][subCol] = this[raw][col]

                subCol++
            }
        }
        subRaw++
    }
    return subMatrix
}

fun MutableList<MutableList<Double>>.getIndex(raw: Int, col: Int): Double = this[raw][col]
fun MutableList<MutableList<Double>>.is2x2(): Boolean = this.size == 2 && this[0].size == 2

// calculate determinant recurve function
fun MutableList<MutableList<Double>>.calculateDeterminant(): Double {
    val determinant: Double
    if (this.is2x2()) {
        determinant = this[0][0] * this[1][1] - this[1][0] * this[0][1]
    } else {
        var cumulativeSum = 0.0
        var sign = 1.0
        for (col in 0 until this[0].size) {
            val subMatrix = this.findSubMatrix(0, col)
            cumulativeSum += this.getIndex(0, col) * sign * subMatrix.calculateDeterminant()
            sign *= -1.0
        }
        determinant = cumulativeSum
    }
    return determinant
}

fun MutableList<MutableList<Double>>.findCoFactor(raw: Int, col: Int): Double {
    return this.findSubMatrix(raw, col).calculateDeterminant() * (-1.0).pow(raw + col)
}

fun MutableList<MutableList<Double>>.adjMatrix(): MutableList<MutableList<Double>> {
    val adjMatrix = MutableList(this.size) { MutableList(this.size) { 0.0 } }
    for (raw in 0 until this.size) {
        for (col in 0 until this[raw].size) {
            adjMatrix[raw][col] = this.findCoFactor(raw, col)
        }
    }
    return adjMatrix
}

fun MutableList<MutableList<Double>>.transpose(): MutableList<MutableList<Double>> {
    val transpose = MutableList(this.size) { MutableList(this.size) { 0.0 } }
    for (raw in 0 until transpose.size) {
        for (col in 0 until transpose[raw].size) {
            transpose[raw][col] = this[col][raw]
        }
    }
    return transpose
}






