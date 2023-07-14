package korlibs.math.geom

// @TODO: Create and expose a half-edge structure here.
data class Triangle(val a: Vector2, val b: Vector2, val c: Vector2) {
    val points = listOf(a, b, c)
}
