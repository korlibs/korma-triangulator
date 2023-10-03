import korlibs.image.color.Colors
import korlibs.korge.Korge
import korlibs.korge.scene.Scene
import korlibs.korge.scene.sceneContainer
import korlibs.korge.view.SContainer
import korlibs.korge.view.align.centerOn
import korlibs.korge.view.graphics
import korlibs.math.geom.Point
import korlibs.math.geom.bezier.Bezier
import korlibs.math.geom.shape.buildVectorPath
import korlibs.math.geom.triangulator.triangulate
import korlibs.math.geom.vector.VectorBuilder
import korlibs.math.geom.vector.VectorPath
import korlibs.math.geom.vector.triangles
import korlibs.math.interpolation.toRatio

suspend fun main() = Korge {
    sceneContainer().changeTo({ MainMyModuleScene() })
}

class MainMyModuleScene : Scene() {
    override suspend fun SContainer.sceneMain() {
        val scale = 4f
        graphics {
            stroke(Colors.WHITE, 1f) {
                triangles(
                    buildVectorPath {
                        circle(Point(0, 0), 50.0 * scale)
                        rect(-20 * scale, -20 * scale, 20 * scale, 20 * scale)
                    }
                    //.processCurves {
                    //    val simpleCurves = it.toSimpleList().map { it.curve }
                    //    //println(simpleCurves.size)
                    //    for (curve in simpleCurves) this.curve(curve)
                    //}
                    .approximate {
                        (it.length / 20f).toInt()
                    }
                    .triangulate()
                )
            }
        }.centerOn(this)
    }
}

fun VectorPath.processCurves(
    process: VectorBuilder.(Bezier) -> Unit
): VectorPath {
    val out = VectorPath()

    visitCmds(
        moveTo = { out.moveTo(it) },
        lineTo = { out.lineTo(it) },
        quadTo = { c, a -> process(out, Bezier(out.lastPos, c, a)) },
        cubicTo = { c1, c2, a -> process(out, Bezier(out.lastPos, c1, c2, a)) },
        close = { out.close() }
    )
    return out
}

fun VectorPath.approximate(
    // @TODO: Provide default. Taking into account the curviness of the curve, and also, taking into account inflection points for cubic
    approximate: (Bezier) -> Int
): VectorPath {
    return processCurves { curve ->
        val out = this
        val count = maxOf(approximate(curve), 3)
        val N = count - 1
        // @TODO: Ensure inflection points are included in the output
        //println(curve.inflections().toList())
        //println("CURVE: $N")
        for (n in 1..N) {
            val ratio = n.toDouble() / N.toDouble()
            out.lineTo(curve[ratio.toRatio()])
        }
    }
}
