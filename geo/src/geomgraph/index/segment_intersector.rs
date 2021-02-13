use super::super::{Edge, Node};
use crate::geomgraph::algorithm::{Intersection, LineIntersector};
use crate::{Coordinate, GeoFloat};

use std::cell::{Ref, RefCell};

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
    fn is_adjacent_segments(&self, i1: usize, i2: usize) -> bool {
        let difference = if i1 > i2 { i1 - i2 } else { i2 - i1 };
        difference == 1
    }

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
    pub fn set_boundary_nodes(
        &mut self,
        boundary_nodes_0: Vec<Node<F>>,
        boundary_nodes_1: Vec<Node<F>>,
    ) {
        // this might be an overzelous assert - JTS doesn't leverage Option types...
        debug_assert!(self.boundary_nodes.is_none());
        self.boundary_nodes = Some([boundary_nodes_0, boundary_nodes_1]);
    }

    pub fn set_is_done_when_proper_intersection(&mut self, new_value: bool) {
        self.is_done_when_proper_intersection = new_value
    }

    pub fn has_proper_intersection(&self) -> bool {
        self.has_proper_intersection
    }

    pub fn has_proper_interior_intersection(&self) -> bool {
        self.has_proper_interior_intersection
    }

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
        // REVIEW: Why is this only applicable when comparing the same edge?
        if *edge0 == *edge1 {
            if self.line_intersector.result() == Intersection::PointIntersection {
                if self.is_adjacent_segments(segment_index_0, segment_index_1) {
                    return true;
                }

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

        false
    }
    pub fn add_intersections(
        &mut self,
        edge0: &RefCell<Edge<F>>,
        segment_index_0: usize,
        edge1: &RefCell<Edge<F>>,
        segment_index_1: usize,
    ) {
        // REVIEW: I *think* we want to compare references here, rather than value equality
        // if *edge0.borrow() == *edge1.borrow() && segment_index_0 == segment_index_1 {
        if edge0.as_ptr() == edge1.as_ptr() && segment_index_0 == segment_index_1 {
            return;
        }

        let p00 = edge0.borrow().coords()[segment_index_0];
        let p01 = edge0.borrow().coords()[segment_index_0 + 1];
        let p10 = edge1.borrow().coords()[segment_index_1];
        let p11 = edge1.borrow().coords()[segment_index_1 + 1];

        self.line_intersector
            .compute_intersection(p00, p01, p10, p11);

        if self.line_intersector.has_intersection() {
            if self.record_isolated {
                edge0.borrow_mut().set_is_isolated(false);
                edge1.borrow_mut().set_is_isolated(false);
            }
            // REVIEW: numIntersections is a private var, never used. Presumably it was for
            //         debugging purposes
            if !self.is_trivial_intersection(
                edge0.borrow(),
                segment_index_0,
                edge1.borrow(),
                segment_index_1,
            ) {
                self.has_intersection = true;

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
                if self.line_intersector.is_proper() {
                    self.proper_intersection_point = Some(self.line_intersector.intersection(0));
                    self.has_proper_intersection = true;

                    if self.is_done_when_proper_intersection {
                        self.is_done = true
                    }

                    if !self.is_boundary_point(&self.line_intersector, &self.boundary_nodes) {
                        self.has_proper_interior_intersection = true
                    }
                }
            }
        }
    }

    fn is_boundary_point(
        &self,
        line_intersector: &Box<dyn LineIntersector<F>>,
        boundary_nodes: &Option<[Vec<Node<F>>; 2]>,
    ) -> bool {
        match boundary_nodes {
            None => false,
            Some(nodes) => {
                self.is_boundary_point_internal(line_intersector, &nodes[0])
                    || self.is_boundary_point_internal(line_intersector, &nodes[1])
            }
        }
    }

    fn is_boundary_point_internal(
        &self,
        line_intersector: &Box<dyn LineIntersector<F>>,
        boundary_nodes: &Vec<Node<F>>,
    ) -> bool {
        boundary_nodes
            .iter()
            .any(|node| line_intersector.is_intersection(node.coordinate()))
    }
}
