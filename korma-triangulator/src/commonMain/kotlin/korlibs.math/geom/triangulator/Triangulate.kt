package korlibs.math.geom.triangulator

import korlibs.math.geom.Point
import korlibs.math.geom.Triangle
import korlibs.math.geom.shape.getPoints2List
import korlibs.math.geom.triangulator.internal.Poly2Tri
import korlibs.math.geom.triangulator.internal.immutable
import korlibs.math.geom.vector.VectorPath

fun triangulate(build: TriangulateBuilder.() -> Unit): List<Triangle> {
    val builder = TriangulateBuilder()
    builder.build()
    return Unit.run { builder.run { finish() } }
}

fun VectorPath.triangulate(): List<Triangle> =
    triangulate { addPolylines(this@triangulate.getPoints2List()) }

class TriangulateBuilder {
    companion object {
        private const val SCALE = 100.0
    }
    private val sweep = Poly2Tri.Sweep()
    private val sc = Poly2Tri.SweepContextExt()

    internal fun Poly2Tri.SweepContextExt.addHoleNew(points: korlibs.math.geom.PointList, closed: Boolean = false) {
        val len = if (closed) points.size else points.size - 1
        for (n in 0 until len) {
            val p1 = points[n]
            val p2 = points[(n + 1) % points.size]
            addEdge(p1 * SCALE, p2 * SCALE)
        }
    }

    fun addPolyline(points: korlibs.math.geom.PointList, closed: Boolean = false) {
        //println("addPolyline: ${points.size}")
        sc.addHoleNew(points, closed = closed)
    }

    fun addPolylines(pointsList: List<korlibs.math.geom.PointList>, closed: Boolean = false) {
        for (points in pointsList) {
            addPolyline(points, closed)
        }
    }

    fun addEdge(a: Point, b: Point) {
        sc.addEdge(a * SCALE, b * SCALE)
    }

    fun Unit.finish(): List<Triangle> {
        sweep.triangulate(sc)
        return sc.getTriangles().map { Triangle(it.a.immutable / SCALE, it.b.immutable / SCALE, it.c.immutable / SCALE) }
    }
}
