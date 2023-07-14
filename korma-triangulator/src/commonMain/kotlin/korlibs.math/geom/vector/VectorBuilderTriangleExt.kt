package korlibs.math.geom.vector

import korlibs.math.geom.Triangle

fun VectorBuilder.triangle(triangle: Triangle) {
    polygon(triangle.points.map { it }, close = true)
}
fun VectorBuilder.triangles(triangles: List<Triangle>) {
    for (triangle in triangles) triangle(triangle)
}
