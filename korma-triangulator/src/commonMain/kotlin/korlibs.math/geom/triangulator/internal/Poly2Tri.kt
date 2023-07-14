package korlibs.math.geom.triangulator.internal

/*
 * Poly2Tri Copyright (c) 2009-2018, Poly2Tri Contributors
 * https://github.com/jhasse/poly2tri
 *
 * Ported by @soywiz from C++ to Kotlin in Jul-2023
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * * Neither the name of Poly2Tri nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

import korlibs.datastructure.linkedHashMapOf
import korlibs.math.geom.Vector2
import kotlin.contracts.ExperimentalContracts
import kotlin.contracts.InvocationKind
import kotlin.contracts.contract
import kotlin.math.absoluteValue

private typealias LinkedList<T> = ArrayList<T>

internal object Poly2Tri {
    const val PI_3DIV4 = 3 * kotlin.math.PI / 4;
    const val PI_DIV2 = 1.57079632679489661923;
    const val EPSILON = 1e-12;

    enum class Orientation { CW, CCW, COLLINEAR }

    val DEBUG = false

    @OptIn(ExperimentalContracts::class)
    inline fun debug(crossinline block: () -> String) {
        contract { callsInPlace(block, InvocationKind.EXACTLY_ONCE) }
        if (DEBUG) println(block())
    }

    /**
     * Formula to calculate signed area
     * Positive if CCW
     * Negative if CW
     * 0 if collinear
     * <pre>
     * A[P1,P2,P3]  =  (x1*y2 - y1*x2) + (x2*y3 - y2*x3) + (x3*y1 - y3*x1)
     *              =  (x1-x3)*(y2-y3) - (y1-y3)*(x2-x3)
     * </pre>
     */
    fun orient2d(pa: Point, pb: Point, pc: Point): Orientation {
        val detleft = (pa.x - pc.x) * (pb.y - pc.y)
        val detright = (pa.y - pc.y) * (pb.x - pc.x)
        val `val` = detleft - detright

        return when {
            //kotlin.math.abs(`val`) < EPSILON -> Orientation.COLLINEAR
            `val` == 0.0 -> Orientation.COLLINEAR
            `val` > 0 -> Orientation.CCW
            else -> Orientation.CW
        }
    }

    fun inScanArea(pa: Point, pb: Point, pc: Point, pd: Point): Boolean {
        val oadb = (pa.x - pb.x)*(pd.y - pb.y) - (pd.x - pb.x)*(pa.y - pb.y)
        if (oadb >= -EPSILON) {
            return false
        }

        val oadc = (pa.x - pc.x)*(pd.y - pc.y) - (pd.x - pc.x)*(pa.y - pc.y)
        if (oadc <= EPSILON) {
            return false
        }
        return true
    }

    // Represents a simple polygon's edge
    class Edge(p1: Point, p2: Point) {
        var p: Point
        var q: Point

        override fun toString(): String = "($p, $q)"

        init {
            when {
                p1.y > p2.y -> {
                    q = p1
                    p = p2
                }
                p1.y == p2.y -> {
                    if (p1.x > p2.x) {
                        q = p1
                        p = p2
                    } else if (p1.x == p2.x) {
                        // Repeat points
                        throw IllegalArgumentException("Edge::Edge: p1 == p2 : $p1 == $p2")
                    } else {
                        q = p2
                        p = p1
                    }
                }
                else -> {
                    q = p2
                    p = p1
                }
            }

            q.edgeList.add(this)
        }
    }

    data class Point(var x: Double = 0.0, var y: Double = 0.0) {
        //override fun toString(): String = "(${x.niceStr(4)}, ${y.niceStr(4)})"
        override fun toString(): String = "($x, $y)"
        // Assuming Edge is a class defined elsewhere
        val edgeList = mutableListOf<Edge>()

        fun setZero() {
            x = 0.0
            y = 0.0
        }

        fun set(x_: Double, y_: Double) {
            x = x_
            y = y_
        }

        operator fun unaryMinus(): Point = Point(-x, -y)

        operator fun plusAssign(p: Point) {
            x += p.x
            y += p.y
        }

        operator fun minusAssign(p: Point) {
            x -= p.x
            y -= p.y
        }

        operator fun timesAssign(a: Double) {
            x *= a
            y *= a
        }

        fun length(): Double = kotlin.math.sqrt(x * x + y * y)

        fun normalize(): Double {
            val len = length()
            x /= len
            y /= len
            return len
        }
    }

    class Triangle(a: Point, b: Point, c: Point) {
        override fun toString(): String = "[$a, $b, $c]"

        val points = arrayOf(a, b, c)
        val a get() = points[0]
        val b get() = points[1]
        val c get() = points[2]
        val neighbors = arrayOfNulls<Triangle?>(3)
        val constrainedEdge = BooleanArray(3) { false }
        val delaunayEdge = BooleanArray(3) { false }
        var interior = false

        fun getPoint(index: Int) = points[index]
        fun getNeighbor(index: Int) = neighbors[index]

        // The point counter-clockwise to given point
        fun pointCW(point: Point): Point {
            return when (point) {
                points[0] -> points[2]
                points[1] -> points[0]
                points[2] -> points[1]
                else -> throw IllegalArgumentException("Invalid point")
            }
        }

        fun pointCCW(point: Point): Point {
            return when (point) {
                points[0] -> points[1]
                points[1] -> points[2]
                points[2] -> points[0]
                else -> throw IllegalArgumentException("Invalid point")
            }
        }

        fun oppositePoint(t: Triangle, p: Point): Point {
            val cw = t.pointCW(p)
            return pointCW(cw!!)!!
        }

        fun markNeighbor(p1: Point, p2: Point, t: Triangle) {
            debug { "MarkNeighbor.p2.t: $p1, $p2, $t" }

            when {
                (p1 == points[2] && p2 == points[1]) || (p1 == points[1] && p2 == points[2]) -> neighbors[0] = t
                (p1 == points[0] && p2 == points[2]) || (p1 == points[2] && p2 == points[0]) -> neighbors[1] = t
                (p1 == points[0] && p2 == points[1]) || (p1 == points[1] && p2 == points[0]) -> neighbors[2] = t
                else -> throw IllegalArgumentException("Points do not match any two points in the Triangle.")
            }
        }

        fun markNeighbor(t: Triangle) {
            debug { "MarkNeighbor.t: $t" }

            when {
                t.contains(points[1], points[2]) -> {
                    neighbors[0] = t
                    t.markNeighbor(points[1], points[2], this)
                }
                t.contains(points[0], points[2]) -> {
                    neighbors[1] = t
                    t.markNeighbor(points[0], points[2], this)
                }
                t.contains(points[0], points[1]) -> {
                    neighbors[2] = t
                    t.markNeighbor(points[0], points[1], this)
                }
            }
        }

        fun markConstrainedEdge(index: Int) {
            constrainedEdge[index] = true;
        }

        fun markConstrainedEdge(edge: Edge) {
            markConstrainedEdge(edge.p, edge.q);
        }

        // Mark edge as constrained
        fun markConstrainedEdge(p: Point, q: Point) {
            when {
                (q == points[0] && p == points[1]) || (q == points[1] && p == points[0]) -> constrainedEdge[2] = true
                (q == points[0] && p == points[2]) || (q == points[2] && p == points[0]) -> constrainedEdge[1] = true
                (q == points[1] && p == points[2]) || (q == points[2] && p == points[1]) -> constrainedEdge[0] = true
            }
        }

        fun index(p: Point): Int {
            return when (p) {
                points[0] -> 0
                points[1] -> 1
                points[2] -> 2
                else -> throw IllegalArgumentException("The point is not part of the triangle.")
            }
        }

        fun edgeIndex(p1: Point, p2: Point): Int {
            return when {
                points[0] == p1 -> {
                    when (p2) {
                        points[1] -> 2
                        points[2] -> 1
                        else -> -1
                    }
                }
                points[1] == p1 -> {
                    when (p2) {
                        points[2] -> 0
                        points[0] -> 2
                        else -> -1
                    }
                }
                points[2] == p1 -> {
                    when (p2) {
                        points[0] -> 1
                        points[1] -> 0
                        else -> -1
                    }
                }
                else -> -1
            }
        }

        // The neighbor across to given point
        fun neighborAcross(point: Point): Triangle? {
            return when (point) {
                points[0] -> neighbors[0]
                points[1] -> neighbors[1]
                else -> neighbors[2]
            }
        }

        // The neighbor clockwise to given point
        fun neighborCW(point: Point): Triangle? {
            return when (point) {
                points[0] -> neighbors[1]
                points[1] -> neighbors[2]
                else -> neighbors[0]
            }
        }

        // The neighbor counter-clockwise to given point
        fun neighborCCW(point: Point): Triangle? {
            return when (point) {
                points[0] -> neighbors[2]
                points[1] -> neighbors[0]
                else -> neighbors[1]
            }
        }

        fun getConstrainedEdgeCCW(p: Point): Boolean {
            return when (p) {
                points[0] -> constrainedEdge[2]
                points[1] -> constrainedEdge[0]
                else -> constrainedEdge[1]
            }
        }

        fun getConstrainedEdgeCW(p: Point): Boolean {
            return when (p) {
                points[0] -> constrainedEdge[1]
                points[1] -> constrainedEdge[2]
                else -> constrainedEdge[0]
            }
        }

        fun setConstrainedEdgeCCW(p: Point, ce: Boolean) {
            when (p) {
                points[0] -> constrainedEdge[2] = ce
                points[1] -> constrainedEdge[0] = ce
                else -> constrainedEdge[1] = ce
            }
        }

        fun setConstrainedEdgeCW(p: Point, ce: Boolean) {
            when (p) {
                points[0] -> constrainedEdge[1] = ce
                points[1] -> constrainedEdge[2] = ce
                else -> constrainedEdge[0] = ce
            }
        }

        fun getDelunayEdgeCCW(p: Point): Boolean {
            return when (p) {
                points[0] -> delaunayEdge[2]
                points[1] -> delaunayEdge[0]
                else -> delaunayEdge[1]
            }
        }

        fun getDelunayEdgeCW(p: Point): Boolean {
            return when (p) {
                points[0] -> delaunayEdge[1]
                points[1] -> delaunayEdge[2]
                else -> delaunayEdge[0]
            }
        }

        fun setDelunayEdgeCCW(p: Point, e: Boolean) {
            when (p) {
                points[0] -> delaunayEdge[2] = e
                points[1] -> delaunayEdge[0] = e
                else -> delaunayEdge[1] = e
            }
        }

        fun setDelunayEdgeCW(p: Point, e: Boolean) {
            when (p) {
                points[0] -> delaunayEdge[1] = e
                points[1] -> delaunayEdge[2] = e
                else -> delaunayEdge[0] = e
            }
        }

        fun contains(p: Point): Boolean = p in points
        fun contains(e: Edge): Boolean = contains(e.p) && contains(e.q)
        fun contains(p: Point, q: Point): Boolean = contains(p) && contains(q)

        // Legalized triangle by rotating clockwise around point(0)
        fun legalize(point: Point) {
            points[1] = points[0]
            points[0] = points[2]
            points[2] = point
        }

        fun legalize(oPoint: Point, nPoint: Point) {
            when {
                oPoint === points[0] -> {
                    points[1] = points[0]
                    points[0] = points[2]
                    points[2] = nPoint
                }
                oPoint === points[1] -> {
                    points[2] = points[1]
                    points[1] = points[0]
                    points[0] = nPoint
                }
                oPoint === points[2] -> {
                    points[0] = points[2]
                    points[2] = points[1]
                    points[1] = nPoint
                }
                else -> throw IllegalArgumentException("The oPoint is not part of the triangle.")
            }
        }

        /**
         * Clears all references to all other triangles and points
         */
        fun clear() {
            for (neighbor in neighbors) {
                neighbor?.clearNeighbor(this)
            }
            clearNeighbors()
            for (i in 0 until points.size) {
                //points[i] = null
            }
        }


        fun clearNeighbor(triangle: Triangle) {
            when (triangle) {
                neighbors[0] -> neighbors[0] = null
                neighbors[1] -> neighbors[1] = null
                else -> neighbors[2] = null
            }
        }

        fun clearNeighbors() {
            neighbors[0] = null
            neighbors[1] = null
            neighbors[2] = null
        }

        fun clearDelunayEdges() {
            delaunayEdge[0] = false
            delaunayEdge[1] = false
            delaunayEdge[2] = false
        }

        fun debugPrint() {
            debug { "${points[0]} ${points[1]} ${points[2]}" }
        }

        fun circumcicleContains(point: Point): Boolean {
            require(isCounterClockwise()) { "Points are not in counter clockwise order" }

            val dx = points[0].x - point.x
            val dy = points[0].y - point.y
            val ex = points[1].x - point.x
            val ey = points[1].y - point.y
            val fx = points[2].x - point.x
            val fy = points[2].y - point.y

            val ap = dx * dx + dy * dy
            val bp = ex * ex + ey * ey
            val cp = fx * fx + fy * fy

            return (dx * (fy * bp - cp * ey) - dy * (fx * bp - cp * ex) + ap * (fx * ey - fy * ex)) < 0
        }

        fun isCounterClockwise(): Boolean {
            return (points[1].x - points[0].x) * (points[2].y - points[0].y) -
                (points[2].x - points[0].x) * (points[1].y - points[0].y) > 0
        }

        companion object {
            fun isDelaunay(triangles: List<Triangle>): Boolean {
                for (triangle in triangles) {
                    for (other in triangles) {
                        if (triangle == other) continue
                        for (i in 0 until 3) {
                            if (triangle.circumcicleContains(other.points[i])) {
                                return false
                            }
                        }
                    }
                }
                return true
            }
        }
    }

    fun icmp(a: Point, b: Point): Int {
        val isZero = a.y.absoluteValue == 0.0 && b.y.absoluteValue == 0.0
        if (!isZero) {
            val res = a.y.compareTo(b.y)
            if (res != 0) return res
        }
        return a.x.compareTo(b.x)
    }

    fun cmp(a: Point, b: Point): Boolean = when {
        a.y < b.y -> true
        a.y == b.y && a.x < b.x -> true
        else -> false
    }

    operator fun Point.plus(b: Point): Point = Point(this.x + b.x, this.y + b.y)
    operator fun Point.minus(b: Point): Point = Point(this.x - b.x, this.y - b.y)
    operator fun Double.times(a: Point): Point = Point(this * a.x, this * a.y)
    fun dot(a: Point, b: Point): Double = a.x * b.x + a.y * b.y
    fun cross(a: Point, b: Point): Double = a.x * b.y - a.y * b.x
    fun cross(a: Point, s: Double): Point = Point(s * a.y, -s * a.x)
    fun cross(s: Double, a: Point): Point = Point(-s * a.y, s * a.x)

    data class Node(
        var point: Point,
        var triangle: Triangle?,
        var next: Node?,
        var prev: Node?,
        var value: Double
    ) {
        constructor(p: Point) : this(p, null, null, null, p.x)
        constructor(p: Point, t: Triangle) : this(p, t, null, null, p.x)
    }


    class AdvancingFront(
        var head: Node,
        var tail: Node,
    ) {
        var searchNode: Node = head

        fun locateNode(x: Double): Node? {
            var node = searchNode

            if (x < node.value) {
                while (node.prev != null) {
                    node = node.prev!!
                    if (x >= node.value) {
                        searchNode = node
                        return node
                    }
                }
            } else {
                while (node.next != null) {
                    node = node.next!!
                    if (x < node.value) {
                        searchNode = node.prev!!
                        return node.prev
                    }
                }
            }
            return null
        }

        fun findSearchNode(x: Double): Node?
        {
            // TODO: implement BST index
            return searchNode;
        }

        fun locatePoint(point: Point): Node? {
            val px = point.x
            var node = findSearchNode(px)
            val nx = node?.point?.x

            if (px == nx) {
                if (point != node?.point) {
                    // We might have two nodes with same x value for a short time
                    if (point == node?.prev?.point) {
                        node = node?.prev
                    } else if (point == node?.next?.point) {
                        node = node?.next
                    } else {
                        error("Invalid")
                    }
                }
            } else if (px < nx!!) {
                while (node?.prev != null) {
                    node = node?.prev
                    if (point == node?.point) {
                        break
                    }
                }
            } else {
                while (node?.next != null) {
                    node = node?.next
                    if (point == node?.point)
                        break
                }
            }
            searchNode = node!!
            return node
        }
    }

    class CDT(private var polylines: List<List<Point>>) {

        private var sweepContext: SweepContext? = null
        private var sweep: Sweep? = null

        init {
            sweepContext = SweepContext(polylines)
            sweep = Sweep()
        }

        fun addHole(polyline: ArrayList<Point>, closed: Boolean) {
            sweepContext?.addHole(polyline, closed)
        }

        fun addPoint(point: Point) {
            sweepContext?.addPoint(point)
        }

        fun triangulate() {
            sweep?.triangulate(sweepContext!!)
        }

        fun getTriangles(): ArrayList<Triangle> {
            return sweepContext?.getTriangles() ?: ArrayList()
        }

        fun getMap(): LinkedList<Triangle> {
            return sweepContext?.getMap() ?: LinkedList()
        }
    }

    open class SweepContextExt(polylines: List<List<Point>> = emptyList()) : SweepContext(polylines) {
        val pointsToIndex = linkedHashMapOf<Vector2, Int>()

        fun addPointToListAndGetIndex(p: Vector2): Int {
            return pointsToIndex.getOrPut(p) {
                points.add(Point(p.xD, p.yD))
                points.size - 1
            }
        }

        fun addEdge(a: Vector2, b: Vector2) {
            if (a == b) return
            val i1 = addPointToListAndGetIndex(a)
            val i2 = addPointToListAndGetIndex(b)
            val p1 = points[i1]
            val p2 = points[i2]
            val edge = Edge(p1, p2)
            this.edgeList.add(edge)
            val index = if (edge.q === p1) i1 else i2
            edgeList.add(edge)
            //println("EDGE[$index]: $p1, $p2  :  $i1, $i2  :  $points_")
        }
    }

    open class SweepContext(polylines: List<List<Point>> = emptyList()) {
        companion object {
            const val K_ALPHA = 0.3
        }

        var edgeList = ArrayList<Edge>()
        var basin = Basin()
        var edgeEvent = EdgeEvent()

        private var triangles = ArrayList<Triangle>()
        private var map = LinkedList<Triangle>()
        protected var points = ArrayList<Point>()

        private var front: AdvancingFront? = null
        private var head: Point? = null
        private var tail: Point? = null

        private var afHead: Node? = null
        private var afMiddle: Node? = null
        private var afTail: Node? = null

        init {
            if (polylines.isNotEmpty()) {
                addHoles(polylines)
            }
        }

        fun addHoles(polylines: List<List<Point>>, closed: Boolean = true) {
            for (polyline in polylines) {
                addHole(polyline, closed)
            }
        }

        fun addHole(polyline: List<Point>, closed: Boolean = true) {
            initEdges(polyline, closed)
            points.addAll(polyline)
        }

        fun addPoint(point: Point) {
            points.add(point)
        }

        fun getTriangles(): ArrayList<Triangle> {
            return triangles
        }

        fun getMap(): LinkedList<Triangle> {
            return map
        }

        internal fun initTriangulation() {
            var xmax = points[0].x
            var xmin = points[0].x
            var ymax = points[0].y
            var ymin = points[0].y

            // Calculate bounds.
            for (point in points) {
                val p = point
                if (p.x > xmax)
                    xmax = p.x
                if (p.x < xmin)
                    xmin = p.x
                if (p.y > ymax)
                    ymax = p.y
                if (p.y < ymin)
                    ymin = p.y
            }

            val dx = K_ALPHA * (xmax - xmin)
            val dy = K_ALPHA * (ymax - ymin)
            head = Point(xmin - dx, ymin - dy)
            tail = Point(xmax + dx, ymin - dy)

            // Sort points along y-axis
            points.sortWith { a, b -> icmp(a, b) }

            for (n in 0 until points.size) {
                debug { "SORTED POINT[$n]: ${points[n]}" }
            }
        }

        private fun initEdges(polyline: List<Point>, closed: Boolean) {
            val numPoints = polyline.size
            for (i in 0 until (if (closed) numPoints else (numPoints - 1))) {
                val j = if (i < numPoints - 1) i + 1 else 0
                edgeList.add(Edge(polyline[i], polyline[j]))
            }
        }

        fun getPoint(index: Int): Point {
            return points[index]
        }

        fun addToMap(triangle: Triangle) {
            map.add(triangle)
        }

        fun locateNode(point: Point): Node? {
            // TODO implement search tree
            return front?.locateNode(point.x)
        }

        fun setHead(p1: Point) {
            head = p1
        }

        fun getHead(): Point? {
            return head
        }

        fun setTail(p1: Point) {
            tail = p1
        }

        fun getTail(): Point? {
            return tail
        }

        fun getPointCount(): Int {
            return points.size
        }

        fun getFront(): AdvancingFront? {
            return front
        }

        fun createAdvancingFront() {

            // Initial triangle
            val triangle = Triangle(points[0], head!!, tail!!)

            map.add(triangle)

            afHead = Node(triangle.getPoint(1), triangle)
            afMiddle = Node(triangle.getPoint(0), triangle)
            afTail = Node(triangle.getPoint(2))
            front = AdvancingFront(afHead!!, afTail!!)

            // TODO: More intuitive if head is middles next and not previous?
            //       so swap head and tail
            afHead!!.next = afMiddle
            afMiddle!!.next = afTail
            afMiddle!!.prev = afHead
            afTail!!.prev = afMiddle
        }

        fun removeNode(node: Node) {
            // Memory is managed automatically in Kotlin, so no equivalent required.
        }

        fun mapTriangleToNodes(t: Triangle) {
            for (i in 0..2) {
                if (t.getNeighbor(i) == null) {
                    val n = front?.locatePoint(t.pointCW(t.getPoint(i)!!))
                    if (n != null)
                        n.triangle = t
                }
            }
        }

        fun removeFromMap(triangle: Triangle) {
            map.remove(triangle)
        }

        fun meshClean(triangle: Triangle) {
            val triangles = mutableListOf(triangle)

            while (triangles.isNotEmpty()){
                val t = triangles.last()
                triangles.removeAt(triangles.size - 1)

                if (!t.interior) {
                    t.interior = true
                    this.triangles.add(t)
                    for (i in 0 until 3) {
                        if (!t.constrainedEdge[i])
                            t.getNeighbor(i)?.let { triangles.add(it) }
                    }
                }
            }
        }

        fun cleanup() {
            // Kotlin has automatic garbage collection, so no need to delete objects

            // but if you need to clear the collections:
            map.clear()
            edgeList.clear()
            triangles.clear()
            // and set the objects to null if needed:
            head = null
            tail = null
            front = null
            afHead = null
            afMiddle = null
            afTail = null
        }
    }

    class Basin {
        var leftNode: Node? = null
        var bottomNode: Node? = null
        var rightNode: Node? = null
        var width: Double = 0.0
        var isLeftHighest: Boolean = false

        fun clear() {
            leftNode = null
            bottomNode = null
            rightNode = null
            width = 0.0
            isLeftHighest = false
        }
    }

    class EdgeEvent {
        var constrainedEdge: Edge? = null
        var isRight: Boolean = false
    }

    open class Sweep {

        private val nodes: MutableList<Node> = ArrayList()

        /**
         * Triangulate
         *
         * @param tcx
         */
        // Triangulate simple polygon with holes
        fun triangulate(tcx: SweepContext) {
            tcx.initTriangulation()
            tcx.createAdvancingFront()

            // Sweep points; build mesh
            sweepPoints(tcx)

            // Clean up
            finalizationPolygon(tcx)
        }

        /**
         * Start sweeping the Y-sorted point set from bottom to top
         *
         * @param tcx
         */
        fun sweepPoints(tcx: SweepContext) {
            for (i in 1 until tcx.getPointCount()) {
                val point = tcx.getPoint(i)
                debug { "SweepPoints.point[$i]: $point" }
                val node = pointEvent(tcx, point)
                for (j in point.edgeList) {
                    debug { "SweepPoints.edge[$i]: $j" }
                    edgeEvent(tcx, j, node)
                }
            }
        }

        /**
         * Find closes node to the left of the new point and
         * create a new triangle. If needed new holes and basins
         * will be filled to.
         *
         * @param tcx
         * @param point
         * @return
         */
        fun pointEvent(tcx: SweepContext, point: Point): Node {
            val nodePtr: Node? = tcx.locateNode(point)
            nodePtr?.point ?: throw RuntimeException("PointEvent - null node")
            nodePtr.next?.point ?: throw RuntimeException("PointEvent - null node")

            val node = nodePtr
            val newNode = newFrontTriangle(tcx, point, node)

            // Only need to check +epsilon since point never have smaller
            // x value than node due to how we fetch nodes from the front
            if (point.x <= node.point.x + EPSILON) {
                fill(tcx, node)
            }

            //tcx.addNode(newNode)

            fillAdvancingFront(tcx, newNode)
            return newNode
        }

        private fun edgeEvent(tcx: SweepContext, edge: Edge, node: Node) {
            tcx.edgeEvent.constrainedEdge = edge
            tcx.edgeEvent.isRight = edge.p.x > edge.q.x

            if (isEdgeSideOfTriangle(node.triangle!!, edge.p, edge.q)) {
                return
            }

            // For now we will do all needed filling
            // TODO: integrate with flip process might give some better performance
            //       but for now this avoid the issue with cases that needs both flips and fills
            fillEdgeEvent(tcx, edge, node)
            edgeEvent(tcx, edge.p, edge.q, node.triangle, edge.q)
        }

        private fun edgeEvent(tcx: SweepContext, ep: Point, eq: Point, triangle: Triangle?, point: Point) {
            debug { "tcx, ep=$ep, eq=$eq, point=$point, triangle=$triangle" }
            var triangle = triangle
            if (triangle == null) {
                throw RuntimeException("EdgeEvent - null triangle")
            }

            debug { " --> [0]" }

            if (isEdgeSideOfTriangle(triangle, ep, eq)) {
                return
            }

            val p1: Point = triangle.pointCCW(point)
            val o1: Orientation = orient2d(eq, p1, ep)

            debug { " --> [1] : $eq, $p1, $ep : $o1" }

            if (o1 == Orientation.COLLINEAR) {
                if (triangle.contains(eq, p1)) {
                    triangle.markConstrainedEdge(eq, p1)
                    // We are modifying the constraint maybe it would be better to
                    // not change the given constraint and just keep a variable for the new constraint
                    tcx.edgeEvent.constrainedEdge!!.q = p1
                    triangle = triangle.neighborAcross(point)
                    edgeEvent(tcx, ep, p1, triangle, p1)
                } else {
                    throw RuntimeException("EdgeEvent - collinear points not supported")
                }
                return
            }

            debug { " --> [2]" }

            val p2: Point = triangle.pointCW(point)
            val o2: Orientation = orient2d(eq, p2, ep)
            if (o2 == Orientation.COLLINEAR) {
                if (triangle.contains(eq, p2)) {
                    triangle.markConstrainedEdge(eq, p2)
                    // We are modifying the constraint maybe it would be better to
                    // not change the given constraint and just keep a variable for the new constraint
                    tcx.edgeEvent.constrainedEdge!!.q = p2
                    triangle = triangle.neighborAcross(point)
                    edgeEvent(tcx, ep, p2, triangle, p2)
                } else {
                    throw RuntimeException("EdgeEvent - collinear points not supported")
                }
                return
            }

            debug { " --> [3]" }

            if (o1 == o2) {
                // Need to decide if we are rotating CW or CCW to get to a triangle
                // that will cross edge
                val otriangle = triangle

                debug { "triangle: ${otriangle.getNeighbor(0)}, ${otriangle.getNeighbor(1)}, ${otriangle.getNeighbor(2)}" }

                triangle = if (o1 == Orientation.CW) triangle.neighborCCW(point) else triangle.neighborCW(point)
                check(triangle != null) {
                    "Triangle is null $otriangle"
                }
                edgeEvent(tcx, ep, eq, triangle, point)
            } else {
                // This triangle crosses constraint so lets flippin start!
                check(triangle != null)
                flipEdgeEvent(tcx, ep, eq, triangle, point)
            }
        }

        private fun newFrontTriangle(tcx: SweepContext, point: Point, node: Node): Node {
            val triangle = Triangle(point, node.point, node.next!!.point)

            triangle.markNeighbor(node.triangle!!)
            tcx.addToMap(triangle)

            val newNode = Node(point)
            nodes.add(newNode)

            newNode.next = node.next
            newNode.prev = node
            node.next!!.prev = newNode
            node.next = newNode

            if (!legalize(tcx, triangle)) {
                tcx.mapTriangleToNodes(triangle)
            }

            return newNode
        }

        private fun fill(tcx: SweepContext, node: Node) {
            val triangle = Triangle(node.prev!!.point, node.point, node.next!!.point)

            // TODO: should copy the constrained_edge value from neighbor triangles
            //       for now constrained_edge values are copied during the legalize
            triangle.markNeighbor(node.prev!!.triangle!!)
            triangle.markNeighbor(node.triangle!!)

            tcx.addToMap(triangle)

            // Update the advancing front
            node.prev!!.next = node.next
            node.next!!.prev = node.prev

            // If it was legalized the triangle has already been mapped
            if (!legalize(tcx, triangle)) {
                tcx.mapTriangleToNodes(triangle)
            }
        }

        private fun legalize(tcx: SweepContext, t: Triangle): Boolean {
            debug { "Legalize: $t" }
            // To legalize a triangle we start by finding if any of the three edges
            // violate the Delaunay condition
            for (i in 0 until 3) {
                debug { "Legalize[a]: $i; ${t.delaunayEdge[i]}" }

                if (t.delaunayEdge[i]) continue

                val ot = t.getNeighbor(i)

                debug { "Legalize[b]: $i; $ot" }


                if (ot == null) continue
                val p = t.getPoint(i)
                val op = ot.oppositePoint(t, p)
                val oi = ot.index(op)

                // If this is a Constrained Edge or a Delaunay Edge(only during recursive legalization)
                // then we should not try to legalize
                if (ot.constrainedEdge[oi] || ot.delaunayEdge[oi]) {
                    t.constrainedEdge[i] = ot.constrainedEdge[oi]
                    continue
                }

                val inside = inCircle(p, t.pointCCW(p), t.pointCW(p), op)

                if (inside) {
                    // Let's mark this shared edge as Delaunay
                    t.delaunayEdge[i] = true
                    ot.delaunayEdge[oi] = true

                    // Let's rotate shared edge one vertex CW to legalize it
                    rotateTrianglePair(t, p, ot, op)

                    // We now got one valid Delaunay Edge shared by two triangles
                    // This gives us 4 new edges to check for Delaunay

                    // Make sure that triangle to node mapping is done only one time for a specific triangle
                    var notLegalized = !legalize(tcx, t)
                    if (notLegalized) {
                        tcx.mapTriangleToNodes(t)
                    }

                    notLegalized = !legalize(tcx, ot)
                    if (notLegalized) {
                        tcx.mapTriangleToNodes(ot)
                    }

                    // Reset the Delaunay edges, since they only are valid Delaunay edges
                    // until we add a new triangle or point.
                    // XXX: need to think about this. Can these edges be tried after we
                    //      return to previous recursive level?
                    t.delaunayEdge[i] = false
                    ot.delaunayEdge[oi] = false

                    // If triangle have been legalized no need to check the other edges since
                    // the recursive legalization will handle those so we can end here.
                    return true
                }

            }
            return false
        }

        private fun inCircle(pa: Point, pb: Point, pc: Point, pd: Point): Boolean {
            val adx = pa.x - pd.x
            val ady = pa.y - pd.y
            val bdx = pb.x - pd.x
            val bdy = pb.y - pd.y

            val adxbdy = adx * bdy
            val bdxady = bdx * ady
            val oabd = adxbdy - bdxady

            if (oabd <= 0) return false

            val cdx = pc.x - pd.x
            val cdy = pc.y - pd.y

            val cdxady = cdx * ady
            val adxcdy = adx * cdy
            val ocad = cdxady - adxcdy

            if (ocad <= 0) return false

            val bdxcdy = bdx * cdy
            val cdxbdy = cdx * bdy

            val alift = adx * adx + ady * ady
            val blift = bdx * bdx + bdy * bdy
            val clift = cdx * cdx + cdy * cdy

            val det = alift * (bdxcdy - cdxbdy) + blift * ocad + clift * oabd

            return det > 0
        }

        fun rotateTrianglePair(t: Triangle, p: Point, ot: Triangle, op: Point) {
            var n1: Triangle? = t.neighborCCW(p)
            var n2: Triangle? = t.neighborCW(p)
            var n3: Triangle? = ot.neighborCCW(op)
            var n4: Triangle? = ot.neighborCW(op)

            val ce1: Boolean = t.getConstrainedEdgeCCW(p)
            val ce2: Boolean = t.getConstrainedEdgeCW(p)
            val ce3: Boolean = ot.getConstrainedEdgeCCW(op)
            val ce4: Boolean = ot.getConstrainedEdgeCW(op)

            val de1: Boolean = t.getDelunayEdgeCCW(p)
            val de2: Boolean = t.getDelunayEdgeCW(p)
            val de3: Boolean = ot.getDelunayEdgeCCW(op)
            val de4: Boolean = ot.getDelunayEdgeCW(op)

            t.legalize(p, op)
            ot.legalize(op, p)

            // Remap delaunay_edge
            ot.setDelunayEdgeCCW(p, de1)
            t.setDelunayEdgeCW(p, de2)
            t.setDelunayEdgeCCW(op, de3)
            ot.setDelunayEdgeCW(op, de4)

            // Remap constrained_edge
            ot.setConstrainedEdgeCCW(p, ce1)
            t.setConstrainedEdgeCW(p, ce2)
            t.setConstrainedEdgeCCW(op, ce3)
            ot.setConstrainedEdgeCW(op, ce4)

            // Remap neighbors
            t.clearNeighbors()
            ot.clearNeighbors()
            n1?.let { ot.markNeighbor(it) }
            n2?.let { t.markNeighbor(it) }
            n3?.let { t.markNeighbor(it) }
            n4?.let { ot.markNeighbor(it) }
            t.markNeighbor(ot)
        }

        private fun fillAdvancingFront(tcx: SweepContext, n: Node) {

            // Fill right holes
            var node = n.next

            while (node != null && node.next != null) {
                // if HoleAngle exceeds 90 degrees then break.
                if (largeHoleDontFill(node)) break
                fill(tcx, node)
                node = node.next
            }

            // Fill left holes
            node = n.prev

            while (node != null && node.prev != null) {
                // if HoleAngle exceeds 90 degrees then break.
                if (largeHoleDontFill(node)) break
                fill(tcx, node)
                node = node.prev
            }

            // Fill right basins
            if (n.next != null && n.next?.next != null) {
                val angle = basinAngle(n)
                if (angle < Poly2Tri.PI_3DIV4) {
                    fillBasin(tcx, n)
                }
            }
        }

        // True if HoleAngle exceeds 90 degrees.
        // LargeHole_DontFill checks if the advancing front has a large hole.
        // A "Large hole" is a triangle formed by a sequence of points in the advancing
        // front where three neighbor points form a triangle.
        // And angle between left-top, bottom, and right-top points is more than 90 degrees.
        // The first part of the algorithm reviews only three neighbor points, e.g. named A, B, C.
        // Additional part of this logic reviews a sequence of 5 points -
        // additionally reviews one point before and one after the sequence of three (A, B, C),
        // e.g. named X and Y.
        // In this case, angles are XBC and ABY and this if angles are negative or more
        // than 90 degrees LargeHole_DontFill returns true.
        // But there is a configuration when ABC has a negative angle but XBC or ABY is less
        // than 90 degrees and positive.
        // Then function LargeHole_DontFill return false and initiates filling.
        // This filling creates a triangle ABC and adds it to the advancing front.
        // But in the case when angle ABC is negative this triangle goes inside the advancing front
        // and can intersect previously created triangles.
        // This triangle leads to making wrong advancing front and problems in triangulation in the future.
        // Looks like such a triangle should not be created.
        // The simplest way to check and fix it is to check an angle ABC.
        // If it is negative LargeHole_DontFill should return true and
        // not initiate creating the ABC triangle in the advancing front.
        // X______A         Y
        //        \        /
        //         \      /
        //          \ B  /
        //           |  /
        //           | /
        //           |/
        //           C
        fun largeHoleDontFill(node: Node): Boolean {

            val nextNode = node.next
            val prevNode = node.prev
            if (!angleExceeds90Degrees(node.point, nextNode?.point!!, prevNode?.point!!))
                return false

            if (angleIsNegative(node.point, nextNode?.point!!, prevNode?.point!!))
                return true

            // Check additional points on front.
            val next2Node = nextNode?.next
            // "..Plus.." because only want angles on same side as point being added.
            if ((next2Node != null) && !angleExceedsPlus90DegreesOrIsNegative(node.point, next2Node.point, prevNode?.point!!))
                return false

            val prev2Node = prevNode?.prev
            // "..Plus.." because only want angles on same side as point being added.
            if ((prev2Node != null) && !angleExceedsPlus90DegreesOrIsNegative(node.point, nextNode?.point!!, prev2Node.point))
                return false

            return true
        }

        private fun angleIsNegative(origin: Point, pa: Point, pb: Point): Boolean {
            val angle = angle(origin, pa, pb)
            return angle < 0
        }

        private fun angleExceeds90Degrees(origin: Point, pa: Point, pb: Point): Boolean {
            val angle = angle(origin, pa, pb)
            return (angle > PI_DIV2) || (angle < -PI_DIV2)
        }

        private fun angleExceedsPlus90DegreesOrIsNegative(origin: Point, pa: Point, pb: Point): Boolean {
            val angle = angle(origin, pa, pb)
            return (angle > PI_DIV2) || (angle < 0)
        }

        private fun angle(origin: Point, pa: Point, pb: Point): Double {
            /* Complex plane
             * ab = cosA +i*sinA
             * ab = (ax + ay*i)(bx + by*i) = (ax*bx + ay*by) + i(ax*by-ay*bx)
             * atan2(y,x) computes the principal value of the argument function
             * applied to the complex number x+iy
             * Where x = ax*bx + ay*by
             *       y = ax*by - ay*bx
             */
            val px = origin.x
            val py = origin.y
            val ax = pa.x - px
            val ay = pa.y - py
            val bx = pb.x - px
            val by = pb.y - py
            val x = ax * by - ay * bx
            val y = ax * bx + ay * by
            return kotlin.math.atan2(x, y)
        }

        private fun holeAngle(node: Node): Double {
            /* Complex plane
             * ab = cosA +i*sinA
             * ab = (ax + ay*i)(bx + by*i) = (ax*bx + ay*by) + i(ax*by-ay*bx)
             * atan2(y,x) computes the principal value of the argument function
             * applied to the complex number x+iy
             * Where x = ax*bx + ay*by
             *       y = ax*by - ay*bx
             */
            val ax = node.next!!.point.x - node.point.x
            val ay = node.next!!.point.y - node.point.y
            val bx = node.prev!!.point.x - node.point.x
            val by = node.prev!!.point.y - node.point.y
            return kotlin.math.atan2(ax * by - ay * bx, ax * bx + ay * by)
        }

        private fun basinAngle(node: Node): Double {
            val ax = node.point.x - node.next!!.next!!.point.x
            val ay = node.point.y - node.next!!.next!!.point.y
            return kotlin.math.atan2(ay, ax)
        }

        private fun fillBasin(tcx: SweepContext, node: Node) {
            tcx.basin.leftNode = if (orient2d(node.point, node.next!!.point, node.next!!.next!!.point) == Orientation.CCW)
                node.next!!.next else node.next

            // Find the bottom and right node
            tcx.basin.bottomNode = tcx.basin.leftNode
            while (tcx.basin.bottomNode!!.next != null && tcx.basin.bottomNode!!.point.y >= tcx.basin.bottomNode!!.next!!.point.y) {
                tcx.basin.bottomNode = tcx.basin.bottomNode!!.next
            }
            if (tcx.basin.bottomNode == tcx.basin.leftNode) {
                // No valid basin
                return
            }

            tcx.basin.rightNode = tcx.basin.bottomNode
            while (tcx.basin.rightNode!!.next != null && tcx.basin.rightNode!!.point.y < tcx.basin.rightNode!!.next!!.point.y) {
                tcx.basin.rightNode = tcx.basin.rightNode!!.next
            }
            if (tcx.basin.rightNode == tcx.basin.bottomNode) {
                // No valid basins
                return
            }

            tcx.basin.width = tcx.basin.rightNode!!.point.x - tcx.basin.leftNode!!.point.x
            tcx.basin.isLeftHighest = tcx.basin.leftNode!!.point.y > tcx.basin.rightNode!!.point.y

            fillBasinReq(tcx, tcx.basin.bottomNode!!)
        }

        private fun fillBasinReq(tcx: SweepContext, node: Node) {
            var node = node
            // if shallow stop filling
            if (isShallow(tcx, node)) {
                return
            }

            fill(tcx, node)

            if (node.prev == tcx.basin.leftNode && node.next == tcx.basin.rightNode) {
                return
            } else if (node.prev == tcx.basin.leftNode) {
                val o = orient2d(node.point, node.next!!.point, node.next!!.next!!.point)
                if (o == Orientation.CW) {
                    return
                }
                node = node.next!!
            } else if (node.next == tcx.basin.rightNode) {
                val o = orient2d(node.point, node.prev!!.point, node.prev!!.prev!!.point)
                if (o == Orientation.CCW) {
                    return
                }
                node = node.prev!!
            } else {
                // Continue with the neighbor node with lowest Y value
                node = if (node.prev!!.point.y < node.next!!.point.y) node.prev!! else node.next!!
            }

            fillBasinReq(tcx, node)
        }

        private fun isShallow(tcx: SweepContext, node: Node): Boolean {
            val height = if (tcx.basin.isLeftHighest) {
                tcx.basin.leftNode!!.point.y - node.point.y
            } else {
                tcx.basin.rightNode!!.point.y - node.point.y
            }

            // if shallow stop filling
            return tcx.basin.width > height
        }

        fun isEdgeSideOfTriangle(triangle: Triangle, ep: Point, eq: Point): Boolean {
            val index = triangle.edgeIndex(ep, eq)

            return if (index != -1) {
                triangle.markConstrainedEdge(index)
                val t: Triangle? = triangle.getNeighbor(index)
                t?.markConstrainedEdge(ep, eq)
                true
            } else {
                false
            }
        }

        private fun fillEdgeEvent(tcx: SweepContext, edge: Edge?, node: Node?) {
            if (tcx.edgeEvent.isRight) {
                fillRightAboveEdgeEvent(tcx, edge, node)
            } else {
                fillLeftAboveEdgeEvent(tcx, edge, node!!)
            }
        }

        private fun fillRightAboveEdgeEvent(tcx: SweepContext, edge: Edge?, node: Node?) {
            var currentNode = node
            while (currentNode!!.next!!.point.x < edge!!.p!!.x) {
                // Check if next node is below the edge
                if (orient2d(edge.q!!, currentNode.next!!.point, edge.p!!) == Orientation.CCW) {
                    fillRightBelowEdgeEvent(tcx, edge, currentNode)
                } else {
                    currentNode = currentNode.next
                }
            }
        }

        private fun fillRightBelowEdgeEvent(tcx: SweepContext, edge: Edge?, node: Node) {
            if (node.point.x < edge!!.p!!.x) {
                if (orient2d(node.point, node.next!!.point, node.next!!.next!!.point) == Orientation.CCW) {
                    // Concave
                    fillRightConcaveEdgeEvent(tcx, edge, node)
                } else {
                    // Convex
                    fillRightConvexEdgeEvent(tcx, edge, node)
                    // Retry this one
                    fillRightBelowEdgeEvent(tcx, edge, node)
                }
            }
        }

        private fun fillRightConcaveEdgeEvent(tcx: SweepContext, edge: Edge?, node: Node) {
            fill(tcx, node.next!!)
            if (node.next!!.point != edge!!.p) {
                // Next above or below edge?
                if (orient2d(edge.q!!, node.next!!.point, edge.p!!) == Orientation.CCW) {
                    // Below
                    if (orient2d(node.point, node.next!!.point, node.next!!.next!!.point) == Orientation.CCW) {
                        // Next is concave
                        fillRightConcaveEdgeEvent(tcx, edge, node)
                    } else {
                        // Next is convex
                    }
                }
            }
        }

        private fun fillRightConvexEdgeEvent(tcx: SweepContext, edge: Edge?, node: Node) {
            // Next concave or convex?
            if (orient2d(node.next!!.point, node.next!!.next!!.point, node.next!!.next!!.next!!.point) == Orientation.CCW) {
                // Concave
                fillRightConcaveEdgeEvent(tcx, edge, node.next!!)
            } else {
                // Convex
                // Next above or below edge?
                if (orient2d(edge!!.q!!, node.next!!.next!!.point, edge.p!!) == Orientation.CCW) {
                    // Below
                    fillRightConvexEdgeEvent(tcx, edge, node.next!!)
                } else {
                    // Above
                }
            }
        }

        private fun fillLeftAboveEdgeEvent(tcx: SweepContext, edge: Edge?, node: Node) {
            var mutableNode = node
            while (mutableNode.prev!!.point.x > edge!!.p!!.x) {
                // Check if next node is below the edge
                if (orient2d(edge.q!!, mutableNode.prev!!.point, edge.p!!) == Orientation.CW) {
                    fillLeftBelowEdgeEvent(tcx, edge, mutableNode)
                } else {
                    mutableNode = mutableNode.prev!!
                }
            }
        }

        private fun fillLeftBelowEdgeEvent(tcx: SweepContext, edge: Edge?, node: Node) {
            var mutableNode = node
            if (mutableNode.point.x > edge!!.p!!.x) {
                if (orient2d(mutableNode.point, mutableNode.prev!!.point, mutableNode.prev!!.prev!!.point) == Orientation.CW) {
                    // Concave
                    fillLeftConcaveEdgeEvent(tcx, edge, mutableNode)
                } else {
                    // Convex
                    fillLeftConvexEdgeEvent(tcx, edge, mutableNode)
                    // Retry this one
                    fillLeftBelowEdgeEvent(tcx, edge, mutableNode)
                }
            }
        }

        private fun fillLeftConcaveEdgeEvent(tcx: SweepContext, edge: Edge?, node: Node) {
            fill(tcx, node.prev!!)
            if (node.prev!!.point != edge!!.p) {
                // Next above or below edge?
                if (orient2d(edge.q!!, node.prev!!.point, edge.p!!) == Orientation.CW) {
                    // Below
                    if (orient2d(node.point, node.prev!!.point, node.prev!!.prev!!.point) == Orientation.CW) {
                        // Next is concave
                        fillLeftConcaveEdgeEvent(tcx, edge, node)
                    } else {
                        // Next is convex
                    }
                }
            }
        }

        private fun fillLeftConvexEdgeEvent(tcx: SweepContext, edge: Edge?, node: Node) {
            // Next concave or convex?
            if (orient2d(node.prev!!.point, node.prev!!.prev!!.point, node.prev!!.prev!!.prev!!.point) == Orientation.CW) {
                // Concave
                fillLeftConcaveEdgeEvent(tcx, edge, node.prev!!)
            } else {
                // Convex
                // Next above or below edge?
                if (orient2d(edge!!.q!!, node.prev!!.prev!!.point, edge.p!!) == Orientation.CW) {
                    // Below
                    fillLeftConvexEdgeEvent(tcx, edge, node.prev!!)
                } else {
                    // Above
                }
            }
        }

        private fun flipEdgeEvent(tcx: SweepContext, ep: Point, eq: Point, t: Triangle, p: Point) {
            val ot = t.neighborAcross(p) ?: throw RuntimeException("FlipEdgeEvent - null neighbor across")
            val op = ot.oppositePoint(t, p)

            if (inScanArea(p, t.pointCCW(p), t.pointCW(p), op)) {
                // Lets rotate shared edge one vertex CW
                rotateTrianglePair(t, p, ot, op)
                tcx.mapTriangleToNodes(t)
                tcx.mapTriangleToNodes(ot)

                if (p == eq && op == ep) {
                    if (eq == tcx.edgeEvent.constrainedEdge!!.q && ep == tcx.edgeEvent.constrainedEdge!!.p) {
                        t.markConstrainedEdge(ep, eq)
                        ot.markConstrainedEdge(ep, eq)
                        legalize(tcx, t)
                        legalize(tcx, ot)
                    } else {
                        // XXX: I think one of the triangles should be legalized here?
                    }
                } else {
                    val o = orient2d(eq, op, ep)
                    val nextFlipTriangle = nextFlipTriangle(tcx, o, t, ot, p, op)
                    flipEdgeEvent(tcx, ep, eq, nextFlipTriangle, p)
                }
            } else {
                val newP = nextFlipPoint(ep, eq, ot, op)
                flipScanEdgeEvent(tcx, ep, eq, t, ot, newP)
                edgeEvent(tcx, ep, eq, t, p)
            }
        }

        private fun nextFlipTriangle(tcx: SweepContext, o: Orientation, t: Triangle, ot: Triangle, p: Point, op: Point): Triangle {
            if (o == Orientation.CCW) {
                // ot is not crossing edge after flip
                val edgeIndex = ot.edgeIndex(p, op)
                ot.delaunayEdge[edgeIndex] = true
                legalize(tcx, ot)
                ot.clearDelunayEdges()
                return t
            }

            // t is not crossing edge after flip
            val edgeIndex = t.edgeIndex(p, op)
            t.delaunayEdge[edgeIndex] = true
            legalize(tcx, t)
            t.clearDelunayEdges()
            return ot
        }

        private fun nextFlipPoint(ep: Point, eq: Point, ot: Triangle, op: Point): Point {
            val o2d = orient2d(eq, op, ep)
            return when (o2d) {
                Orientation.CW -> ot.pointCCW(op)
                Orientation.CCW -> ot.pointCW(op)
                else -> throw UnsupportedOperationException("Opposing point on constrained edge")
            }
        }

        private fun flipScanEdgeEvent(tcx: SweepContext, ep: Point, eq: Point, flipTriangle: Triangle,
                                      t: Triangle, p: Point) {
            val ot = t.neighborAcross(p) ?: throw RuntimeException("FlipScanEdgeEvent - null neighbor across")
            val op = ot.oppositePoint(t, p) ?: throw RuntimeException("FlipScanEdgeEvent - null opposing point")

            val p1 = flipTriangle.pointCCW(eq)
            val p2 = flipTriangle.pointCW(eq)
            if (p1 == null || p2 == null) {
                throw RuntimeException("FlipScanEdgeEvent - null on either of points")
            }

            if (inScanArea(eq, p1, p2, op)) {
                // flip with new edge op->eq
                flipEdgeEvent(tcx, eq, op, ot, op)
                // TODO: Actually I just figured out that it should be possible to
                //       improve this by getting the next ot and op before the above
                //       flip and continue the flipScanEdgeEvent here
                // set new ot and op here and loop back to inScanArea test
                // also need to set a new flipTriangle first
                // Turns out at first glance that this is somewhat complicated
                // so it will have to wait.
            } else {
                val newP = nextFlipPoint(ep, eq, ot, op)
                flipScanEdgeEvent(tcx, ep, eq, flipTriangle, ot, newP)
            }
        }

        fun finalizationPolygon(tcx: SweepContext) {
            // Get an Internal triangle to start with
            var t: Triangle? = tcx.getFront()!!.head.next?.triangle
            val p: Point? = tcx.getFront()!!.head.next?.point
            while (t != null && !t.getConstrainedEdgeCW(p!!)) {
                t = t.neighborCCW(p)
            }

            // Collect interior triangles constrained by edges
            t?.let { tcx.meshClean(it) }
        }
    }
}