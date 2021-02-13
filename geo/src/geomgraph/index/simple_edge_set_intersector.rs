use super::super::Edge;
use super::{EdgeSetIntersector, SegmentIntersector};
use crate::GeoFloat;

use std::cell::RefCell;
use std::rc::Rc;

pub(crate) struct SimpleEdgeSetIntersector {
    overlap_count: u32,
}

impl<F> EdgeSetIntersector<F> for SimpleEdgeSetIntersector
where
    F: GeoFloat,
{
    fn compute_intersections(
        &mut self,
        edges: &[Rc<RefCell<Edge<F>>>],
        segment_intersector: &mut SegmentIntersector<F>,
        test_all_segments: bool,
    ) {
        self.overlap_count = 0;

        for edge0 in edges.iter() {
            for edge1 in edges.iter() {
                // TODO: I expect this to explode if test_all_segments is true, due to the RefCell
                //       being borrowed twice simultaneously
                //
                if test_all_segments || edge0 != edge1 {
                    self.compute_intersects(edge0, edge1, segment_intersector);
                }
            }
        }
    }

    fn compute_intersections_testing_all_segments(
        &mut self,
        edges0: &[Rc<RefCell<Edge<F>>>],
        edges1: &[Rc<RefCell<Edge<F>>>],
        segment_intersector: &mut SegmentIntersector<F>,
    ) {
        self.overlap_count = 0;

        for edge0 in edges0 {
            for edge1 in edges1 {
                self.compute_intersects(edge0, edge1, segment_intersector);
            }
        }
    }
}

impl SimpleEdgeSetIntersector {
    pub fn new() -> Self {
        SimpleEdgeSetIntersector { overlap_count: 0 }
    }

    fn compute_intersects<'a, F: GeoFloat>(
        &mut self,
        edge0: &Rc<RefCell<Edge<F>>>,
        edge1: &Rc<RefCell<Edge<F>>>,
        segment_intersector: &mut SegmentIntersector<F>,
    ) {
        let edge0_coords_len = edge0.borrow().coords().len() - 1;
        let edge1_coords_len = edge1.borrow().coords().len() - 1;
        for i0 in 0..edge0_coords_len {
            for i1 in 0..edge1_coords_len {
                segment_intersector.add_intersections(edge0, i0, edge1, i1);
            }
        }
    }
}
