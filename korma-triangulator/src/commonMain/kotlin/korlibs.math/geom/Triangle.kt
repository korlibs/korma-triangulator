package korlibs.math.geom

// @TODO: Create and expose a half-edge structure here.
data class Triangle(val a: Vector2D, val b: Vector2D, val c: Vector2D) {
    val points = listOf(a, b, c)
}
