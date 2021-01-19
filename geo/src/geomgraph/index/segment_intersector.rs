use super::super::{Edge, Node};
use crate::geomgraph::algorithm::{Intersection, LineIntersector};
use crate::{Coordinate, GeoFloat};

use std::cell::{Ref, RefCell};

// JTS: /**
// JTS:  * Computes the intersection of line segments,
// JTS:  * and adds the intersection to the edges containing the segments.
// JTS:  *
// JTS:  * @version 1.7
// JTS:  */
// JTS: public class SegmentIntersector
// JTS: {
/// Computes the intersection of line segments,
/// and adds the intersection to the edges containing the segments.
pub(crate) struct SegmentIntersector<F>
where
    F: GeoFloat,
{
    // TODO is it worth making this generic?
    // Though JTS leaves this abstract - we might consider hard coding it to a RobustLineIntersector
    // TODO benchmark to see if there is an appreciable perf difference
    line_intersector: Box<dyn LineIntersector<F>>,
    record_isolated: bool,
    has_intersection: bool,
    include_proper: bool,
    // CLEANUP: can we get rid of has_proper_intersection and just do proper_intersection_point.is_some()?
    proper_intersection_point: Option<Coordinate<F>>,
    has_proper_intersection: bool,
    has_proper_interior_intersection: bool,
    is_done_when_proper_intersection: bool,
    is_done: bool,
    boundary_nodes: Option<[Vec<Node<F>>; 2]>,
}

impl<F> SegmentIntersector<F>
where
    F: GeoFloat,
{
    // JTS:
    // JTS:   public static boolean isAdjacentSegments(int i1, int i2)
    // JTS:   {
    // JTS:     return Math.abs(i1 - i2) == 1;
    // JTS:   }
    fn is_adjacent_segments(&self, i1: usize, i2: usize) -> bool {
        let difference = if i1 > i2 { i1 - i2 } else { i2 - i1 };
        difference == 1
    }

    // JTS:
    // JTS:   /**
    // JTS:    * These variables keep track of what types of intersections were
    // JTS:    * found during ALL edges that have been intersected.
    // JTS:    */
    // JTS:   private boolean hasIntersection = false;
    // JTS:   private boolean hasProper = false;
    // JTS:   private boolean hasProperInterior = false;
    // JTS:   // the proper intersection point found
    // JTS:   private Coordinate properIntersectionPoint = null;
    // JTS:
    // JTS:   private LineIntersector li;
    // JTS:   private boolean includeProper;
    // JTS:   private boolean recordIsolated;
    // JTS:   private boolean isSelfIntersection;
    // JTS:   //private boolean intersectionFound;
    // JTS:   private int numIntersections = 0;
    // JTS:
    // JTS:   // testing only
    // JTS:   public int numTests = 0;
    // JTS:
    // JTS:   private Collection[] bdyNodes;
    // JTS:   private boolean isDone = false;
    // JTS:   private boolean isDoneWhenProperInt = false;
    // JTS:
    // JTS:
    // JTS:   public SegmentIntersector(LineIntersector li,  boolean includeProper, boolean recordIsolated)
    // JTS:   {
    // JTS:     this.li = li;
    // JTS:     this.includeProper = includeProper;
    // JTS:     this.recordIsolated = recordIsolated;
    // JTS:   }
    pub fn new(
        line_intersector: Box<dyn LineIntersector<F>>,
        include_proper: bool,
        record_isolated: bool,
    ) -> SegmentIntersector<F> {
        SegmentIntersector {
            line_intersector,
            include_proper,
            record_isolated,
            has_intersection: false,
            has_proper_intersection: false,
            has_proper_interior_intersection: false,
            proper_intersection_point: None,
            is_done: false,
            is_done_when_proper_intersection: false,
            boundary_nodes: None,
        }
    }
    // JTS:
    // JTS:   public void setBoundaryNodes( Collection bdyNodes0,
    // JTS:                               Collection bdyNodes1)
    // JTS:   {
    // JTS:       bdyNodes = new Collection[2];
    // JTS:       bdyNodes[0] = bdyNodes0;
    // JTS:       bdyNodes[1] = bdyNodes1;
    // JTS:   }
    pub fn set_boundary_nodes(
        &mut self,
        boundary_nodes_0: Vec<Node<F>>,
        boundary_nodes_1: Vec<Node<F>>,
    ) {
        // this might be an overzelous assert - JTS doesn't leverage Option types...
        debug_assert!(self.boundary_nodes.is_none());
        self.boundary_nodes = Some([boundary_nodes_0, boundary_nodes_1]);
    }

    // JTS:   public void setIsDoneIfProperInt(boolean isDoneWhenProperInt) {
    // JTS: 	  this.isDoneWhenProperInt = isDoneWhenProperInt;
    // JTS:   }
    pub fn set_is_done_when_proper_intersection(&mut self, new_value: bool) {
        self.is_done_when_proper_intersection = new_value
    }

    // JTS:
    // JTS:   public boolean isDone() {
    // JTS: 	  return isDone;
    // JTS:   }
    // JTS:   /**
    // JTS:    * @return the proper intersection point, or <code>null</code> if none was found
    // JTS:    */
    // JTS:   public Coordinate getProperIntersectionPoint()  {    return properIntersectionPoint;  }
    // JTS:
    // JTS:   public boolean hasIntersection() { return hasIntersection; }
    // JTS:   /**
    // JTS:    * A proper intersection is an intersection which is interior to at least two
    // JTS:    * line segments.  Note that a proper intersection is not necessarily
    // JTS:    * in the interior of the entire Geometry, since another edge may have
    // JTS:    * an endpoint equal to the intersection, which according to SFS semantics
    // JTS:    * can result in the point being on the Boundary of the Geometry.
    // JTS:    */
    // JTS:   public boolean hasProperIntersection() { return hasProper; }
    pub fn has_proper_intersection(&self) -> bool {
        self.has_proper_intersection
    }

    // JTS:   /**
    // JTS:    * A proper interior intersection is a proper intersection which is <b>not</b>
    // JTS:    * contained in the set of boundary nodes set for this SegmentIntersector.
    // JTS:    */
    // JTS:   public boolean hasProperInteriorIntersection() { return hasProperInterior; }
    pub fn has_proper_interior_intersection(&self) -> bool {
        self.has_proper_interior_intersection
    }

    // JTS:
    // JTS:
    // JTS:   /**
    // JTS:    * A trivial intersection is an apparent self-intersection which in fact
    // JTS:    * is simply the point shared by adjacent line segments.
    // JTS:    * Note that closed edges require a special check for the point shared by the beginning
    // JTS:    * and end segments.
    // JTS:    */
    // JTS:   private boolean isTrivialIntersection(Edge e0, int segIndex0, Edge e1, int segIndex1)
    // JTS:   {
    /// A trivial intersection is an apparent self-intersection which in fact is simply the point
    /// shared by adjacent line segments.  Note that closed edges require a special check for the
    /// point shared by the beginning and end segments.
    fn is_trivial_intersection(
        &self,
        edge0: Ref<'_, Edge<F>>,
        segment_index_0: usize,
        edge1: Ref<'_, Edge<F>>,
        segment_index_1: usize,
    ) -> bool {
        // JTS:     if (e0 == e1) {
        // JTS:       if (li.getIntersectionNum() == 1) {
        // REVIEW: Why is this only applicable when comparing the same edge?
        if *edge0 == *edge1 {
            if self.line_intersector.result() == Intersection::PointIntersection {
                // JTS:         if (isAdjacentSegments(segIndex0, segIndex1))
                // JTS:           return true;
                if self.is_adjacent_segments(segment_index_0, segment_index_1) {
                    return true;
                }

                // JTS:         if (e0.isClosed()) {
                // JTS:           int maxSegIndex = e0.getNumPoints() - 1;
                // JTS:           if (    (segIndex0 == 0 && segIndex1 == maxSegIndex)
                // JTS:               ||  (segIndex1 == 0 && segIndex0 == maxSegIndex) ) {
                // JTS:             return true;
                // JTS:           }
                // JTS:         }
                if edge0.is_closed() {
                    let max_segment_index = edge0.coords().len() - 1;
                    if (segment_index_0 == 0 && segment_index_1 == max_segment_index)
                        || (segment_index_1 == 0 && segment_index_0 == max_segment_index)
                    {
                        return true;
                    }
                }
            }
        }
        // JTS:       }
        // JTS:     }

        // JTS:     return false;
        // JTS:   }
        false
    }
    // JTS:
    // JTS:   /**
    // JTS:    * This method is called by clients of the EdgeIntersector class to test for and add
    // JTS:    * intersections for two segments of the edges being intersected.
    // JTS:    * Note that clients (such as MonotoneChainEdges) may choose not to intersect
    // JTS:    * certain pairs of segments for efficiency reasons.
    // JTS:    */
    // JTS:   public void addIntersections(
    // JTS:     Edge e0,  int segIndex0,
    // JTS:     Edge e1,  int segIndex1
    // JTS:      )
    // JTS:   {
    pub fn add_intersections(
        &mut self,
        edge0: &RefCell<Edge<F>>,
        segment_index_0: usize,
        edge1: &RefCell<Edge<F>>,
        segment_index_1: usize,
    ) {
        // JTS:     if (e0 == e1 && segIndex0 == segIndex1) return;
        // REVIEW: I *think* we want to compare references here, rather than value equality
        // if *edge0.borrow() == *edge1.borrow() && segment_index_0 == segment_index_1 {
        if edge0.as_ptr() == edge1.as_ptr() && segment_index_0 == segment_index_1 {
            return;
        }

        // JTS: numTests++;
        // JTS:     Coordinate p00 = e0.getCoordinates()[segIndex0];
        // JTS:     Coordinate p01 = e0.getCoordinates()[segIndex0 + 1];
        // JTS:     Coordinate p10 = e1.getCoordinates()[segIndex1];
        // JTS:     Coordinate p11 = e1.getCoordinates()[segIndex1 + 1];
        let p00 = edge0.borrow().coords()[segment_index_0];
        let p01 = edge0.borrow().coords()[segment_index_0 + 1];
        let p10 = edge1.borrow().coords()[segment_index_1];
        let p11 = edge1.borrow().coords()[segment_index_1 + 1];

        // JTS:     li.computeIntersection(p00, p01, p10, p11);
        self.line_intersector
            .compute_intersection(p00, p01, p10, p11);

        // JTS: //if (li.hasIntersection() && li.isProper()) Debug.println(li);
        // JTS:     /**
        // JTS:      *  Always record any non-proper intersections.
        // JTS:      *  If includeProper is true, record any proper intersections as well.
        // JTS:      */
        // JTS:     if (li.hasIntersection()) {
        if self.line_intersector.has_intersection() {
            // JTS:       if (recordIsolated) {
            // JTS:         e0.setIsolated(false);
            // JTS:         e1.setIsolated(false);
            // JTS:       }
            if self.record_isolated {
                edge0.borrow_mut().set_is_isolated(false);
                edge1.borrow_mut().set_is_isolated(false);
            }
            // JTS:       //intersectionFound = true;
            // REVIEW: numIntersections is a private var, never used. Presumably it was for
            //         debugging purposes
            // JTS:       numIntersections++;
            // JTS:       // if the segments are adjacent they have at least one trivial intersection,
            // JTS:       // the shared endpoint.  Don't bother adding it if it is the
            // JTS:       // only intersection.
            // JTS:       if (! isTrivialIntersection(e0, segIndex0, e1, segIndex1)) {
            if !self.is_trivial_intersection(
                edge0.borrow(),
                segment_index_0,
                edge1.borrow(),
                segment_index_1,
            ) {
                // JTS:         hasIntersection = true;
                self.has_intersection = true;

                // JTS:         if (includeProper || ! li.isProper() ) {
                // JTS: //Debug.println(li);
                // JTS:           e0.addIntersections(li, segIndex0, 0);
                // JTS:           e1.addIntersections(li, segIndex1, 1);
                // JTS:         }
                if self.include_proper || !self.line_intersector.is_proper() {
                    edge0.borrow_mut().add_intersections(
                        &self.line_intersector,
                        segment_index_0,
                        0,
                    );
                    edge1.borrow_mut().add_intersections(
                        &self.line_intersector,
                        segment_index_1,
                        1,
                    );
                }
                // JTS:         if (li.isProper()) {
                if self.line_intersector.is_proper() {
                    // JTS:           properIntersectionPoint = li.getIntersection(0).copy();
                    // JTS:           hasProper = true;
                    self.proper_intersection_point = Some(self.line_intersector.intersection(0));
                    self.has_proper_intersection = true;

                    // JTS:           if (isDoneWhenProperInt) {
                    // JTS:         	  isDone = true;
                    // JTS:           }
                    if self.is_done_when_proper_intersection {
                        self.is_done = true
                    }

                    // JTS:           if (! isBoundaryPoint(li, bdyNodes))
                    // JTS:             hasProperInterior = true;
                    if !self.is_boundary_point(&self.line_intersector, &self.boundary_nodes) {
                        self.has_proper_interior_intersection = true
                    }
                    // JTS:         }
                    // JTS:         //if (li.isCollinear())
                    // JTS:           //hasCollinear = true;
                    // JTS:       }
                }
            }
            // JTS:     }
            // JTS:   }
        }
    }

    // JTS:   private boolean isBoundaryPoint(LineIntersector li, Collection[] bdyNodes)
    // JTS:   {
    fn is_boundary_point(
        &self,
        line_intersector: &Box<dyn LineIntersector<F>>,
        boundary_nodes: &Option<[Vec<Node<F>>; 2]>,
    ) -> bool {
        // JTS:     if (bdyNodes == null) return false;
        // JTS:     if (isBoundaryPointInternal(li, bdyNodes[0])) return true;
        // JTS:     if (isBoundaryPointInternal(li, bdyNodes[1])) return true;
        // JTS:     return false;
        match boundary_nodes {
            None => false,
            Some(nodes) => {
                self.is_boundary_point_internal(line_intersector, &nodes[0])
                    || self.is_boundary_point_internal(line_intersector, &nodes[1])
            }
        }
        // JTS:   }
    }

    // JTS:   private boolean isBoundaryPointInternal(LineIntersector li, Collection bdyNodes)
    // JTS:   {
    fn is_boundary_point_internal(
        &self,
        line_intersector: &Box<dyn LineIntersector<F>>,
        boundary_nodes: &Vec<Node<F>>,
    ) -> bool {
        // JTS:     for (Iterator i = bdyNodes.iterator(); i.hasNext(); ) {
        // JTS:       Node node = (Node) i.next();
        // JTS:       Coordinate pt = node.getCoordinate();
        // JTS:       if (li.isIntersection(pt)) return true;
        // JTS:     }
        // JTS:     return false;
        boundary_nodes
            .iter()
            .any(|node| line_intersector.is_intersection(node.coordinate()))
        // JTS:   }
    }
    // JTS: }
}
