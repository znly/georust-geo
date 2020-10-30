use super::super::{Edge, Float};
use super::SegmentIntersector;

use std::cell::RefCell;

pub(crate) trait EdgeSetIntersector<F>
where
    F: Float,
{
    // JTS: /**
    // JTS:  * Computes all self-intersections between edges in a set of edges,
    // JTS:  * allowing client to choose whether self-intersections are computed.
    // JTS:  *
    // JTS:  * @param edges a list of edges to test for intersections
    // JTS:  * @param si the SegmentIntersector to use
    // JTS:  * @param testAllSegments true if self-intersections are to be tested as well
    // JTS:  */
    // JTS: abstract public void computeIntersections(List edges, SegmentIntersector si, boolean testAllSegments);
    fn compute_intersections(
        &mut self,
        edges: &[RefCell<Edge<F>>],
        segment_intersector: &mut SegmentIntersector<F>,
        test_all_segments: bool,
    );

    // JTS: /**
    // JTS:   * Computes all mutual intersections between two sets of edges.
    // JTS:   */
    // JTS:  abstract public void computeIntersections(List edges0, List edges1, SegmentIntersector si);
    fn compute_intersections_testing_all_segments(
        &mut self,
        edges0: &[RefCell<Edge<F>>],
        edges1: &[RefCell<Edge<F>>],
        segment_intersector: &mut SegmentIntersector<F>,
    );
}
