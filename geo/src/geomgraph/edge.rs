use super::algorithm::LineIntersector;
use super::{Dimensions, EdgeIntersectionList, Label, Position};
use crate::{Coordinate, GeoFloat};

// REVIEW: kind of weird circular dependency on relate module from geomgraph
use crate::algorithm::relate::IntersectionMatrix;

// TODO: investigate how isEqual should be implented - not sure it makes sense
// to derive equality, since it compares a bunch of vec's
#[derive(Debug)]
pub(crate) struct Edge<F: GeoFloat> {
    coords: Vec<Coordinate<F>>,
    is_isolated: bool,
    edge_intersections: EdgeIntersectionList<F>,
    label: Label,
}

/// Graph Component methods
impl<F: GeoFloat> Edge<F> {
    pub(crate) fn label(&self) -> &Label {
        &self.label
    }

    pub(crate) fn label_mut(&mut self) -> &mut Label {
        &mut self.label
    }
}

impl<F: GeoFloat> Edge<F> {
    pub fn update_intersection_matrix(label: &Label, intersection_matrix: &mut IntersectionMatrix) {
        intersection_matrix.set_at_least_if_valid(
            label.location(0, Position::On),
            label.location(1, Position::On),
            Dimensions::OneDimensional,
        );

        if label.is_area() {
            intersection_matrix.set_at_least_if_valid(
                label.location(0, Position::Left),
                label.location(1, Position::Left),
                Dimensions::TwoDimensional,
            );
            intersection_matrix.set_at_least_if_valid(
                label.location(0, Position::Right),
                label.location(1, Position::Right),
                Dimensions::TwoDimensional,
            );
        }
    }
}

impl<F: GeoFloat> Edge<F> {
    pub fn coords(&self) -> &[Coordinate<F>] {
        &self.coords
    }

    pub fn is_isolated(&self) -> bool {
        self.is_isolated
    }
    pub fn set_is_isolated(&mut self, new_value: bool) {
        self.is_isolated = new_value;
    }

    pub fn new(coords: Vec<Coordinate<F>>, label: Label) -> Edge<F> {
        Edge {
            coords,
            label: label,
            is_isolated: true,
            edge_intersections: EdgeIntersectionList::new(),
        }
    }
    pub fn coordinate(&self) -> Option<&Coordinate<F>> {
        self.coords.get(0)
    }

    pub fn edge_intersections(&self) -> &EdgeIntersectionList<F> {
        &self.edge_intersections
    }

    pub fn edge_intersections_mut(&mut self) -> &mut EdgeIntersectionList<F> {
        &mut self.edge_intersections
    }

    pub fn add_edge_intersection_list_endpoints(&mut self) {
        let max_segment_index = self.coords().len() - 1;
        let first_coord = self.coords()[0];
        let max_coord = self.coords()[max_segment_index];
        self.edge_intersections_mut().add(first_coord, 0, F::zero());
        self.edge_intersections_mut()
            .add(max_coord, max_segment_index, F::zero());
    }

    pub fn is_closed(&self) -> bool {
        self.coords().first() == self.coords().last()
    }

    /// Adds EdgeIntersections for one or both intersections found for a segment of an edge to the
    /// edge intersection list.
    pub fn add_intersections(
        &mut self,
        line_intersector: &Box<dyn LineIntersector<F>>,
        segment_index: usize,
        geom_index: usize,
    ) {
        for i in 0..line_intersector.intersection_num() {
            self.add_intersection(line_intersector, segment_index, geom_index, i);
        }
    }

    /// Add an EdgeIntersection for intersection intIndex.
    /// An intersection that falls exactly on a vertex of the edge is normalized to use the higher
    /// of the two possible segmentIndexes
    pub fn add_intersection(
        &mut self,
        line_intersector: &Box<dyn LineIntersector<F>>,
        segment_index: usize,
        geom_index: usize,
        intersection_index: usize,
    ) {
        let intersection_coord = line_intersector.intersection(intersection_index);

        let mut normalized_segment_index = segment_index;
        let mut distance = line_intersector.edge_distance(geom_index, intersection_index);

        let next_segment_index = normalized_segment_index + 1;

        if next_segment_index < self.coords.len() {
            let next_coord = self.coords[next_segment_index];
            if intersection_coord == next_coord {
                normalized_segment_index = next_segment_index;
                distance = F::zero();
            }
        }
        self.edge_intersections
            .add(intersection_coord, normalized_segment_index, distance);
    }
}

impl<F: GeoFloat> std::cmp::PartialEq for Edge<F> {
    fn eq(&self, other: &Edge<F>) -> bool {
        if self.coords.len() != other.coords.len() {
            return false;
        }

        let mut is_equal_forward = true;
        let mut is_equal_reverse = true;
        let coords_len = self.coords.len();
        for i in 0..coords_len {
            if self.coords()[i] != other.coords()[i] {
                is_equal_forward = false;
            }

            if self.coords()[i] != other.coords()[coords_len - i - 1] {
                is_equal_reverse = false;
            }

            if !is_equal_forward && !is_equal_reverse {
                return false;
            }
        }

        return true;
    }
}

impl<F: GeoFloat> Edge<F> {}
