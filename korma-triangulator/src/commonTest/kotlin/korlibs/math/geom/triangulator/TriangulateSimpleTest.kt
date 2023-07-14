package korlibs.math.geom.triangulator

import korlibs.io.async.suspendTest
import korlibs.math.geom.Point
import korlibs.math.geom.shape.buildVectorPath
import korlibs.math.geom.triangulator.internal.Poly2Tri
import kotlin.test.Test
import kotlin.test.assertEquals

class TriangulateSimpleTest {
    @Test
    fun testTriangulation3() = suspendTest {
        buildVectorPath {
            circle(Point(500, 500), 500f)
            //rect(0, 0, 100, 100)
            line(Point(400, 400), Point(400, 600))
            line(Point(400, 600), Point(600, 400))
        }.triangulate()
    }

    @Test
    fun testTriangulation2() = suspendTest {
        fun TriangulateBuilder.addSquare() {
            addEdge(Point(0, 0), Point(100, 0))
            addEdge(Point(100, 0), Point(100, 100))
            addEdge(Point(100, 100), Point(0, 100))
            addEdge(Point(0, 100), Point(0, 0))
        }

        assertEquals(2, triangulate {
            addSquare()
        }.size)

        assertEquals(6, triangulate {
            addSquare()
            addEdge(Point(40, 40), Point(40, 60))
        }.size)
        assertEquals(8, triangulate {
            addSquare()
            addEdge(Point(40, 40), Point(40, 60))
            addEdge(Point(40,60), Point(60, 40))
        }.size)

        println(triangulate {
            addSquare()
            addEdge(Point(40, 40), Point(40, 60))
            addEdge(Point(40, 60), Point(60, 40))
        })

        //NativeImageContext2d(200, 200) {
        //    stroke(Colors.GREEN, 2f) {
        //        triangles(triangulate {
        //            addSquare()
        //            addEdge(Point(40, 40), Point(40, 60))
        //            addEdge(Point(40, 60), Point(60, 40))
        //        })
        //    }
        //}.writeTo("/tmp/image.png".uniVfs, PNG)
    }

    @Test
    fun testTriangulation() = suspendTest {
        val sweep = Poly2Tri.Sweep()
        val sc = Poly2Tri.SweepContextExt()
        sc.addEdge(Point(0, 0), Point(100, 0))
        sc.addEdge(Point(100, 0), Point(100, 100))
        sc.addEdge(Point(100, 100), Point(0, 100))
        sc.addEdge(Point(0, 100), Point(0, 0))
        sc.addEdge(Point(40, 40), Point(40, 60))
        //sc.addEdge(Point(40, 60), Point(60, 40))

        //sc.addHole(listOf(MPoint(0, 0), MPoint(100, 0), MPoint(100, 100), MPoint(0, 100)), closed = true)
        //sc.addPoint(MPoint(40, 40))
        //sc.addPoint(MPoint(70, 40))
        //sc.addHole(listOf(MPoint(40, 40), MPoint(40, 60)), closed = false)
        //sc.addHole(listOf(MPoint(100, 100), MPoint(0, 100)))

        sweep.triangulate(sc)
        println(sc.getTriangles().toList())
    }
}