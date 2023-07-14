package korlibs.math.geom.triangulator

import korlibs.io.async.*
import korlibs.io.file.*
import korlibs.io.file.std.*
import korlibs.math.geom.*
import korlibs.math.geom.triangulator.internal.Poly2Tri
import kotlin.test.*

class TriangulationTestBedCommonTest {

    @Test
    fun testSimpleAddEdge() {
        val sweep = Poly2Tri.Sweep()
        val context = Poly2Tri.SweepContextExt()
        context.addEdge(Point(0, 0), Point(100, 0))
        context.addEdge(Point(100, 0), Point(100, 100))
        context.addEdge(Point(100, 100), Point(0, 100))
        context.addEdge(Point(0, 100), Point(0, 0))
        context.addEdge(Point(40, 40), Point(40, 60))
        sweep.triangulate(context)
        println(context.getTriangles())
    }

    @Test fun test2() = testFile(resourcesVfs["data/2.dat"], 58)
    @Test fun test_bird() = testFile(resourcesVfs["data/bird.dat"], 273)
    @Test fun test_custom() = testFile(resourcesVfs["data/custom.dat"], 3)
    @Test fun test_deadly_quad() = testFile(resourcesVfs["data/deadly_quad.dat"], 2)
    @Test fun test_deadly_quad_relaxed() = testFile(resourcesVfs["data/deadly_quad_relaxed.dat"], 2)
    @Test fun test_debug() = testFile(resourcesVfs["data/debug.dat"], 198)
    @Test fun test_debug2() = testFile(resourcesVfs["data/debug2.dat"], 9998)
    @Test fun test_diamond() = testFile(resourcesVfs["data/diamond.dat"], 8)
    @Test fun test_dude() = testFile(resourcesVfs["data/dude.dat"], 106)
    @Test fun test_funny() = testFile(resourcesVfs["data/funny.dat"], 98)
    @Test fun test_kzer() = testFile(resourcesVfs["data/kzer-za.dat"], 206)
    @Test fun test_nazca_heron() = testFile(resourcesVfs["data/nazca_heron.dat"], 1034)
    @Test fun test_nazca_monkey() = testFile(resourcesVfs["data/nazca_monkey.dat"], 1202)
    @Test fun test_polygon_test_01() = testFile(resourcesVfs["data/polygon_test_01.dat"], 10)
    @Test fun test_polygon_test_02() = testFile(resourcesVfs["data/polygon_test_02.dat"], 11)
    @Test fun test_polygon_test_03() = testFile(resourcesVfs["data/polygon_test_03.dat"], 5)
    //@Test fun test_sketchup() = testFile(resourcesVfs["data/sketchup.dat"])
    @Test fun test_stalactite() = testFile(resourcesVfs["data/stalactite.dat"], 13)
    @Test fun test_star() = testFile(resourcesVfs["data/star.dat"], 8)
    @Test fun test_steiner() = testFile(resourcesVfs["data/steiner.dat"], 10)
    @Test fun test_strange() = testFile(resourcesVfs["data/strange.dat"], 14)
    @Test fun test_tank() = testFile(resourcesVfs["data/tank.dat"], 53)
    @Test fun test_test() = testFile(resourcesVfs["data/test.dat"], 4)


    fun testFile(file: VfsFile, ntriangles: Int = 0) = suspendTest {
        val pointsList = listPointListFromString(file.readString())
        val sweep = Poly2Tri.Sweep()
        val sweepContext = Poly2Tri.SweepContext()
        for (points in pointsList) {
            //println("points=$points")
            val pp = points.points.map { Poly2Tri.Point(it.x, it.y) }
            if (points.steiner) {
                for (p in pp) sweepContext.addPoint(p)
            } else {
                sweepContext.addHole(pp)
            }
        }
        sweep.triangulate(sweepContext)
        val tris = sweepContext.getTriangles()
        assertEquals(ntriangles, tris.size)
        //println("tris=$tris")
        //val sweep = Poly2Tri.SweepContextExt(points)
        //val tris = sweep.triangulate().getTriangles()
        //for (tri in tris) Unit
    }
}
