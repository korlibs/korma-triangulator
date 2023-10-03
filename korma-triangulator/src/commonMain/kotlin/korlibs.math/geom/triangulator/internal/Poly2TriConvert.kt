package korlibs.math.geom.triangulator.internal

import korlibs.math.geom.Triangle
import korlibs.math.geom.Vector2D

internal val Poly2Tri.Point.immutable: Vector2D get() = Vector2D(x, y)
internal val Poly2Tri.Triangle.immutable: Triangle get() = Triangle(a.immutable, b.immutable, c.immutable)
