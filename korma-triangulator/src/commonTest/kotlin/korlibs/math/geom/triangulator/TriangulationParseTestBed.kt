package korlibs.math.geom.triangulator

import korlibs.math.geom.*

data class Poly2TriPointList(val points: List<MPoint>, val steiner: Boolean)

fun listPointListFromString(string: String): List<Poly2TriPointList> {
    val holes = arrayListOf<Poly2TriPointList>()
    var currentPoints = arrayListOf<MPoint>()
    var steiner = false

    fun split() {
        if (currentPoints.isNotEmpty()) {
            holes += Poly2TriPointList(currentPoints, steiner)
            currentPoints = arrayListOf()
        }
        steiner = false
    }

    for (str in string.split("\n").map { it.trim() }) {
        if (str.isEmpty()) continue
        if (str.contains(" ") || str.contains("\t")) {
            val (p1, p2) = str.trim().split(Regex("\\s+"))
            val cp = MPoint(p1.toDouble(), p2.toDouble())
            if (currentPoints.lastOrNull() != cp) {
                currentPoints.add(cp)
            }
        } else {
            when (str) {
                "HOLE" -> split()
                "STEINER" -> {
                    split()
                    steiner = true
                }
                else -> error("Unsupported '$str'")
            }
        }
    }
    split()

    return holes
}
