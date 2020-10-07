// JTS: import org.locationtech.jts.geom.Coordinate;
// JTS: import org.locationtech.jts.geom.PrecisionModel;
// JTS: import org.locationtech.jts.io.WKTWriter;
// JTS: import org.locationtech.jts.util.Assert;
use geo_types::Coordinate;

#[derive(PartialEq, Eq, Clone, Copy)]
pub enum Intersection {
    /// Indicates that line segments do not intersect
    NoIntersection,

    /// Indicates that line segments intersect in a single point
    PointIntersection,

    /// Indicates that line segments intersect in a line segment
    CollinearIntersection,
}

// JTS: /**
// JTS:  * A <code>LineIntersector</code> is an algorithm that can both test whether
// JTS:  * two line segments intersect and compute the intersection point(s)
// JTS:  * if they do.
// JTS:  * <p>
// JTS:  * There are three possible outcomes when determining whether two line segments intersect:
// JTS:  * <ul>
// JTS:  * <li>{@link #NO_INTERSECTION} - the segments do not intersect
// JTS:  * <li>{@link #POINT_INTERSECTION} - the segments intersect in a single point
// JTS:  * <li>{@link #COLLINEAR_INTERSECTION} - the segments are collinear and they intersect in a line segment
// JTS:  * </ul>
// JTS:  * For segments which intersect in a single point, the point may be either an endpoint
// JTS:  * or in the interior of each segment.
// JTS:  * If the point lies in the interior of both segments,
// JTS:  * this is termed a <i>proper intersection</i>.
// JTS:  * The method {@link #isProper()} test for this situation.
// JTS:  * <p>
// JTS:  * The intersection point(s) may be computed in a precise or non-precise manner.
// JTS:  * Computing an intersection point precisely involves rounding it
// JTS:  * via a supplied {@link PrecisionModel}.
// JTS:  * <p>
// JTS:  * LineIntersectors do not perform an initial envelope intersection test
// JTS:  * to determine if the segments are disjoint.
// JTS:  * This is because this class is likely to be used in a context where
// JTS:  * envelope overlap is already known to occur (or be likely).
// JTS:  *
// JTS:  * @version 1.7
// JTS:  */
// JTS: public abstract class LineIntersector
// JTS: {

pub(crate) trait LineIntersector<F: num_traits::Float> {
    // JTS: /**
    // JTS:  * These are deprecated, due to ambiguous naming
    // JTS:  */
    // JTS:   public final static int DONT_INTERSECT = 0;
    // JTS:   public final static int DO_INTERSECT = 1;
    // JTS:   public final static int COLLINEAR = 2;
    // JTS:
    // JTS:   /**
    // JTS:    * Indicates that line segments do not intersect
    // JTS:    */
    // JTS:   public final static int NO_INTERSECTION = 0;
    // JTS:
    // JTS:   /**
    // JTS:    * Indicates that line segments intersect in a single point
    // JTS:    */
    // JTS:   public final static int POINT_INTERSECTION = 1;
    // JTS:
    // JTS:   /**
    // JTS:    * Indicates that line segments intersect in a line segment
    // JTS:    */
    // JTS:   public final static int COLLINEAR_INTERSECTION = 2;
    // JTS:
    // JTS:   /**
    // JTS:    * Computes the "edge distance" of an intersection point p along a segment.
    // JTS:    * The edge distance is a metric of the point along the edge.
    // JTS:    * The metric used is a robust and easy to compute metric function.
    // JTS:    * It is <b>not</b> equivalent to the usual Euclidean metric.
    // JTS:    * It relies on the fact that either the x or the y ordinates of the
    // JTS:    * points in the edge are unique, depending on whether the edge is longer in
    // JTS:    * the horizontal or vertical direction.
    // JTS:    * <p>
    // JTS:    * NOTE: This function may produce incorrect distances
    // JTS:    *  for inputs where p is not precisely on p1-p2
    // JTS:    * (E.g. p = (139,9) p1 = (139,10), p2 = (280,1) produces distance 0.0, which is incorrect.
    // JTS:    * <p>
    // JTS:    * My hypothesis is that the function is safe to use for points which are the
    // JTS:    * result of <b>rounding</b> points which lie on the line,
    // JTS:    * but not safe to use for <b>truncated</b> points.
    // JTS:    */
    // JTS:   public static double computeEdgeDistance(
    // JTS:         Coordinate p,
    // JTS:         Coordinate p0,
    // JTS:         Coordinate p1)
    // JTS:   {
    fn compute_edge_distance(&self, p: Coordinate<F>, p0: Coordinate<F>, p1: Coordinate<F>) -> F {
        // JTS:     double dx = Math.abs(p1.x - p0.x);
        // JTS:     double dy = Math.abs(p1.y - p0.y);
        let dx = p1.x - p0.x;
        let dy = p1.y - p0.y;

        // JTS:     double dist = -1.0;   // sentinel value
        let mut dist: F;
        // JTS:     if (p.equals(p0)) {
        // JTS:       dist = 0.0;
        // JTS:     }
        if p == p0 {
            dist = F::zero();
        // JTS:     else if (p.equals(p1)) {
        // JTS:       if (dx > dy)
        // JTS:         dist = dx;
        // JTS:       else
        // JTS:         dist = dy;
        // JTS:     }
        } else if p == p1 {
            if dx > dy {
                dist = dx;
            } else {
                dist = dy;
            }
        // JTS:     else {
        // JTS:       double pdx = Math.abs(p.x - p0.x);
        // JTS:       double pdy = Math.abs(p.y - p0.y);
        // JTS:       if (dx > dy)
        // JTS:         dist = pdx;
        // JTS:       else
        // JTS:         dist = pdy;
        // JTS:       // <FIX>
        // JTS:       // hack to ensure that non-endpoints always have a non-zero distance
        // JTS:       if (dist == 0.0 && ! p.equals(p0))
        // JTS:       {
        // JTS:         dist = Math.max(pdx, pdy);
        // JTS:       }
        // JTS:     }
        } else {
            let pdx = (p.x - p0.x).abs();
            let pdy = (p.y - p0.y).abs();
            if dx > dy {
                dist = pdx;
            } else {
                dist = pdy;
            }
            // hack to ensure that non-endpoints always have a non-zero distance
            if dist == F::zero() && p != p0 {
                dist = pdx.max(pdy);
            }
        }
        // JTS:     Assert.isTrue(! (dist == 0.0 && ! p.equals(p0)), "Bad distance calculation");
        // JTS:     return dist;
        assert!(!(dist == F::zero() && p != p0));
        dist
        // JTS:   }
    }
    // JTS:
    // JTS:   /**
    // JTS:    * This function is non-robust, since it may compute the square of large numbers.
    // JTS:    * Currently not sure how to improve this.
    // JTS:    */
    // JTS:   public static double nonRobustComputeEdgeDistance(
    // JTS:         Coordinate p,
    // JTS:         Coordinate p1,
    // JTS:         Coordinate p2)
    // JTS:   {
    // JTS:     double dx = p.x - p1.x;
    // JTS:     double dy = p.y - p1.y;
    // JTS:     double dist = Math.sqrt(dx * dx + dy * dy);   // dummy value
    // JTS:     Assert.isTrue(! (dist == 0.0 && ! p.equals(p1)), "Invalid distance calculation");
    // JTS:     return dist;
    // JTS:   }
    // JTS:

    // JTS:   protected int result;
    fn get_result(&self) -> Intersection;
    fn set_result(&mut self, intersection: Intersection);

    // JTS:   protected Coordinate[][] inputLines = new Coordinate[2][2];
    fn get_input_lines(&self) -> [[Coordinate<F>; 2]; 2];
    fn set_input_lines(&mut self, line: usize, start_or_end: usize, coord: Coordinate<F>);

    // JTS:   protected Coordinate[] intPt = new Coordinate[2];
    // JTS:   /**
    // JTS:    * The indexes of the endpoints of the intersection lines, in order along
    // JTS:    * the corresponding line
    // JTS:    */
    // JTS:   protected int[][] intLineIndex;
    // JTS:   protected boolean isProper;
    // JTS:   protected Coordinate pa;
    // JTS:   protected Coordinate pb;
    // JTS:   /**
    // JTS:    * If makePrecise is true, computed intersection coordinates will be made precise
    // JTS:    * using Coordinate#makePrecise
    // JTS:    */
    // JTS:   protected PrecisionModel precisionModel = null;
    // JTS: //public int numIntersects = 0;
    // JTS:
    // JTS:   public LineIntersector() {
    // JTS:     intPt[0] = new Coordinate();
    // JTS:     intPt[1] = new Coordinate();
    // JTS:     // alias the intersection points for ease of reference
    // JTS:     pa = intPt[0];
    // JTS:     pb = intPt[1];
    // JTS:     result = 0;
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Force computed intersection to be rounded to a given precision model
    // JTS:    * @param precisionModel
    // JTS:    * @deprecated use <code>setPrecisionModel</code> instead
    // JTS:    */
    // JTS:   public void setMakePrecise(PrecisionModel precisionModel)
    // JTS:   {
    // JTS:     this.precisionModel = precisionModel;
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Force computed intersection to be rounded to a given precision model.
    // JTS:    * No getter is provided, because the precision model is not required to be specified.
    // JTS:    * @param precisionModel
    // JTS:    */
    // JTS:   public void setPrecisionModel(PrecisionModel precisionModel)
    // JTS:   {
    // JTS:     this.precisionModel = precisionModel;
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Gets an endpoint of an input segment.
    // JTS:    *
    // JTS:    * @param segmentIndex the index of the input segment (0 or 1)
    // JTS:    * @param ptIndex the index of the endpoint (0 or 1)
    // JTS:    * @return the specified endpoint
    // JTS:    */
    // JTS:   public Coordinate getEndpoint(int segmentIndex, int ptIndex)
    // JTS:   {
    // JTS:     return inputLines[segmentIndex][ptIndex];
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Compute the intersection of a point p and the line p1-p2.
    // JTS:    * This function computes the boolean value of the hasIntersection test.
    // JTS:    * The actual value of the intersection (if there is one)
    // JTS:    * is equal to the value of <code>p</code>.
    // JTS:    */
    // JTS:   public abstract void computeIntersection(
    // JTS:         Coordinate p,
    // JTS:         Coordinate p1, Coordinate p2);
    // JTS:
    // JTS:   protected boolean isCollinear() {
    // JTS:     return result == COLLINEAR_INTERSECTION;
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Computes the intersection of the lines p1-p2 and p3-p4.
    // JTS:    * This function computes both the boolean value of the hasIntersection test
    // JTS:    * and the (approximate) value of the intersection point itself (if there is one).
    // JTS:    */
    // JTS:   public void computeIntersection(
    // JTS:                 Coordinate p1, Coordinate p2,
    // JTS:                 Coordinate p3, Coordinate p4) {
    fn compute_intersection(
        &mut self,
        p1: Coordinate<F>,
        p2: Coordinate<F>,
        p3: Coordinate<F>,
        p4: Coordinate<F>,
    ) {
        // JTS:     inputLines[0][0] = p1;
        // JTS:     inputLines[0][1] = p2;
        // JTS:     inputLines[1][0] = p3;
        // JTS:     inputLines[1][1] = p4;
        // JTS:     result = computeIntersect(p1, p2, p3, p4);
        // JTS: //numIntersects++;
        // JTS:   }
        // TODO: What do we need input_lines for? This API feels overly stateful, but I'm going to
        //       match JTS for now.
        self.set_input_lines(0, 0, p1);
        self.set_input_lines(0, 1, p2);
        self.set_input_lines(1, 0, p3);
        self.set_input_lines(1, 1, p4);

        // TODO: This API feels overly stateful, but I'm going to match JTS for now.
        let result = self.compute_intersect(p1, p2, p3, p4);
        self.set_result(result);
    }

    // JTS:
    // JTS:   protected abstract int computeIntersect(
    // JTS:                 Coordinate p1, Coordinate p2,
    // JTS:                 Coordinate q1, Coordinate q2);
    // JTS:
    fn compute_intersect(
        &mut self,
        p1: Coordinate<F>,
        p2: Coordinate<F>,
        q1: Coordinate<F>,
        q2: Coordinate<F>,
    ) -> Intersection;

    // JTS: /*
    // JTS:   public String toString() {
    // JTS:     String str = inputLines[0][0] + "-"
    // JTS:          + inputLines[0][1] + " "
    // JTS:          + inputLines[1][0] + "-"
    // JTS:          + inputLines[1][1] + " : "
    // JTS:                + getTopologySummary();
    // JTS:     return str;
    // JTS:   }
    // JTS: */
    // JTS:
    // JTS:   public String toString() {
    // JTS:     return WKTWriter.toLineString(inputLines[0][0], inputLines[0][1]) + " - "
    // JTS:     + WKTWriter.toLineString(inputLines[1][0], inputLines[1][1])
    // JTS:                  + getTopologySummary();
    // JTS:   }
    // JTS:
    // JTS:   private String getTopologySummary()
    // JTS:   {
    // JTS:     StringBuilder catBuilder = new StringBuilder();
    // JTS:     if (isEndPoint()) catBuilder.append(" endpoint");
    // JTS:     if (isProper) catBuilder.append(" proper");
    // JTS:     if (isCollinear()) catBuilder.append(" collinear");
    // JTS:     return catBuilder.toString();
    // JTS:   }
    // JTS:
    // JTS:   protected boolean isEndPoint() {
    // JTS:     return hasIntersection() && !isProper;
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Tests whether the input geometries intersect.
    // JTS:    *
    // JTS:    * @return true if the input geometries intersect
    // JTS:    */
    // JTS:   public boolean hasIntersection() {
    // JTS:     return result != NO_INTERSECTION;
    // JTS:   }

    /// Tests whether the input geometries intersect.
    fn has_intersection(&self) -> bool {
        self.get_result() != Intersection::NoIntersection
    }

    // JTS:
    // JTS:   /**
    // JTS:    * Returns the number of intersection points found.  This will be either 0, 1 or 2.
    // JTS:    *
    // JTS:    * @return the number of intersection points found (0, 1, or 2)
    // JTS:    */
    // JTS:   public int getIntersectionNum() { return result; }
    /// Returns the number of intersection points found.  This will be either 0, 1 or 2.
    ///
    /// @return the number of intersection points found (0, 1, or 2)
    fn get_intersection_num(&self) -> usize {
        match self.get_result() {
            Intersection::NoIntersection => 0,
            Intersection::PointIntersection => 1,
            Intersection::CollinearIntersection => 2,
        }
    }

    // JTS:
    // JTS:   /**
    // JTS:    * Returns the intIndex'th intersection point
    // JTS:    *
    // JTS:    * @param intIndex is 0 or 1
    // JTS:    *
    // JTS:    * @return the intIndex'th intersection point
    // JTS:    */
    // JTS:   public Coordinate getIntersection(int intIndex)  { return intPt[intIndex]; }
    fn get_intersection(&self, intersection_index: usize) -> Coordinate<F>;

    // JTS:
    // JTS:   protected void computeIntLineIndex() {
    // JTS:     if (intLineIndex == null) {
    // JTS:       intLineIndex = new int[2][2];
    // JTS:       computeIntLineIndex(0);
    // JTS:       computeIntLineIndex(1);
    // JTS:     }
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Test whether a point is a intersection point of two line segments.
    // JTS:    * Note that if the intersection is a line segment, this method only tests for
    // JTS:    * equality with the endpoints of the intersection segment.
    // JTS:    * It does <b>not</b> return true if
    // JTS:    * the input point is internal to the intersection segment.
    // JTS:    *
    // JTS:    * @return true if the input point is one of the intersection points.
    // JTS:    */
    // JTS:   public boolean isIntersection(Coordinate pt) {
    // JTS:     for (int i = 0; i < result; i++) {
    // JTS:       if (intPt[i].equals2D(pt)) {
    // JTS:         return true;
    // JTS:       }
    // JTS:     }
    // JTS:     return false;
    // JTS:   }
    /// Test whether a point is a intersection point of two line segments.
    ///
    /// Note that if the intersection is a line segment, this method only tests for equality with
    /// the endpoints of the intersection segment.  
    ///
    /// It does <b>not</b> return true if the input point is internal to the intersection segment.
    ///
    /// @return true if the input point is one of the intersection points.
    fn is_intersection(&self, coord: &Coordinate<F>) -> bool {
        for i in 0..self.get_intersection_num() {
            if self.get_intersection(i) == *coord {
                return true;
            }
        }
        return false;
    }

    // JTS:
    // JTS:   /**
    // JTS:    * Tests whether either intersection point is an interior point of one of the input segments.
    // JTS:    *
    // JTS:    * @return <code>true</code> if either intersection point is in the interior of one of the input segments
    // JTS:    */
    // JTS:   public boolean isInteriorIntersection()
    // JTS:   {
    // JTS:     if (isInteriorIntersection(0)) return true;
    // JTS:     if (isInteriorIntersection(1)) return true;
    // JTS:     return false;
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Tests whether either intersection point is an interior point of the specified input segment.
    // JTS:    *
    // JTS:    * @return <code>true</code> if either intersection point is in the interior of the input segment
    // JTS:    */
    // JTS:   public boolean isInteriorIntersection(int inputLineIndex)
    // JTS:   {
    // JTS:     for (int i = 0; i < result; i++) {
    // JTS:       if (! (   intPt[i].equals2D(inputLines[inputLineIndex][0])
    // JTS:              || intPt[i].equals2D(inputLines[inputLineIndex][1]) )) {
    // JTS:         return true;
    // JTS:       }
    // JTS:     }
    // JTS:     return false;
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Tests whether an intersection is proper.
    // JTS:    * <br>
    // JTS:    * The intersection between two line segments is considered proper if
    // JTS:    * they intersect in a single point in the interior of both segments
    // JTS:    * (e.g. the intersection is a single point and is not equal to any of the
    // JTS:    * endpoints).
    // JTS:    * <p>
    // JTS:    * The intersection between a point and a line segment is considered proper
    // JTS:    * if the point lies in the interior of the segment (e.g. is not equal to
    // JTS:    * either of the endpoints).
    // JTS:    *
    // JTS:    * @return true if the intersection is proper
    // JTS:    */
    // JTS:   public boolean isProper() {
    // JTS:     return hasIntersection() && isProper;
    // JTS:   }
    fn is_proper(&self) -> bool;
    // JTS:
    // JTS:   /**
    // JTS:    * Computes the intIndex'th intersection point in the direction of
    // JTS:    * a specified input line segment
    // JTS:    *
    // JTS:    * @param segmentIndex is 0 or 1
    // JTS:    * @param intIndex is 0 or 1
    // JTS:    *
    // JTS:    * @return the intIndex'th intersection point in the direction of the specified input line segment
    // JTS:    */
    // JTS:   public Coordinate getIntersectionAlongSegment(int segmentIndex, int intIndex) {
    // JTS:     // lazily compute int line array
    // JTS:     computeIntLineIndex();
    // JTS:     return intPt[intLineIndex[segmentIndex][intIndex]];
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Computes the index (order) of the intIndex'th intersection point in the direction of
    // JTS:    * a specified input line segment
    // JTS:    *
    // JTS:    * @param segmentIndex is 0 or 1
    // JTS:    * @param intIndex is 0 or 1
    // JTS:    *
    // JTS:    * @return the index of the intersection point along the input segment (0 or 1)
    // JTS:    */
    // JTS:   public int getIndexAlongSegment(int segmentIndex, int intIndex) {
    // JTS:     computeIntLineIndex();
    // JTS:     return intLineIndex[segmentIndex][intIndex];
    // JTS:   }
    // JTS:
    // JTS:   protected void computeIntLineIndex(int segmentIndex) {
    // JTS:     double dist0 = getEdgeDistance(segmentIndex, 0);
    // JTS:     double dist1 = getEdgeDistance(segmentIndex, 1);
    // JTS:     if (dist0 > dist1) {
    // JTS:       intLineIndex[segmentIndex][0] = 0;
    // JTS:       intLineIndex[segmentIndex][1] = 1;
    // JTS:     }
    // JTS:     else {
    // JTS:       intLineIndex[segmentIndex][0] = 1;
    // JTS:       intLineIndex[segmentIndex][1] = 0;
    // JTS:     }
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Computes the "edge distance" of an intersection point along the specified input line segment.
    // JTS:    *
    // JTS:    * @param segmentIndex is 0 or 1
    // JTS:    * @param intIndex is 0 or 1
    // JTS:    *
    // JTS:    * @return the edge distance of the intersection point
    // JTS:    */
    // JTS:   public double getEdgeDistance(int segmentIndex, int intIndex) {
    // JTS:     double dist = computeEdgeDistance(intPt[intIndex], inputLines[segmentIndex][0],
    // JTS:         inputLines[segmentIndex][1]);
    // JTS:     return dist;
    // JTS:   }
    fn get_edge_distance(&self, segment_index: usize, intersection_index: usize) -> F {
        self.compute_edge_distance(
            self.get_intersection(intersection_index),
            self.get_input_lines()[segment_index][0],
            self.get_input_lines()[segment_index][1],
        )
    }
    // JTS: }
}
