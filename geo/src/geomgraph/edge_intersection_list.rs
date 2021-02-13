use super::{Edge, EdgeIntersection};
use crate::{Coordinate, GeoFloat};

use std::collections::BTreeSet;

#[derive(PartialEq, Debug)]
pub(crate) struct EdgeIntersectionList<F>
where
    F: GeoFloat,
{
    node_set: BTreeSet<EdgeIntersection<F>>,
}

impl<F> EdgeIntersectionList<F>
where
    F: GeoFloat,
{
    pub fn new() -> EdgeIntersectionList<F> {
        // REVIEW: JTS holds a circular reference to edge.
        // Let's see if we can just pass it in to the methods that need it instead.
        EdgeIntersectionList {
            node_set: BTreeSet::new(),
        }
    }

    /// Adds an intersection into the list, if it isn't already there.
    /// The input segmentIndex and dist are expected to be normalized.
    /// @return the EdgeIntersection found or added
    pub fn add(&mut self, intersection_point: Coordinate<F>, segment_index: usize, dist: F) {
        let edge_intersection = EdgeIntersection::new(intersection_point, segment_index, dist);
        // BTreeSet only updates the element if it's not alread present
        self.node_set.insert(edge_intersection);

        // Note: the JTS implementation returns the new EdgeIntersection, but it seems unused.
        // Returning it would require some reference gymnastics, so I'm going to omit it until such
        // a time as its needed.
    }
}

impl<'a, F: GeoFloat> IntoIterator for &'a EdgeIntersectionList<F> {
    type Item = &'a EdgeIntersection<F>;
    type IntoIter = std::collections::btree_set::Iter<'a, EdgeIntersection<F>>;

    fn into_iter(self) -> Self::IntoIter {
        self.node_set.iter()
    }
}

impl<F> EdgeIntersectionList<F>
where
    F: GeoFloat,
{
    // REVIEW: Note, we pass in the edge rather than maintain a RefCell to it
    // pub fn add_endpoints(&mut self, edge: &mut Edge<F>) {
    //     let max_segment_index = edge.coords().len() - 1;
    //     self.add(edge.coords()[0], 0, F::zero());
    //     self.add(edge.coords()[max_segment_index], max_segment_index, F::zero());
    // }
}
