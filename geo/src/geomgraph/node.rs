use super::{Dimensions, EdgeEnd, EdgeEndBundleStar, Label, Location};
use crate::{Coordinate, GeoFloat};

// weird circular dependency from GeomGraph to IntersectionMatrix
use crate::algorithm::relate::IntersectionMatrix;

#[derive(Debug, Clone)]
pub(crate) struct Node<F>
where
    F: GeoFloat,
{
    coordinate: Coordinate<F>,
    label: Label,
}

impl<F: GeoFloat> Node<F> {
    pub(crate) fn label(&self) -> &Label {
        &self.label
    }

    pub(crate) fn label_mut(&mut self) -> &mut Label {
        &mut self.label
    }

    pub(crate) fn is_isolated(&self) -> bool {
        self.label.geometry_count() == 1
    }
}

impl<F> Node<F>
where
    F: GeoFloat,
{
    pub fn new(coordinate: Coordinate<F>) -> Node<F> {
        Node {
            coordinate,
            label: Label::empty_line(),
        }
    }

    pub fn coordinate(&self) -> &Coordinate<F> {
        &self.coordinate
    }

    pub fn set_label_on_location(&mut self, geom_index: usize, location: Location) {
        self.label.set_on_location(geom_index, location)
    }

    pub fn set_label_boundary(&mut self, geom_index: usize) {
        let new_location = match self.label.on_location(geom_index) {
            Some(Location::Boundary) => Location::Interior,
            Some(Location::Interior) => Location::Boundary,
            None | Some(Location::Exterior) => Location::Boundary,
        };
        self.label.set_on_location(geom_index, new_location);
    }
}

impl<F: GeoFloat> Node<F> {
    // from JTS#GraphComponent - seems like only node uses this impl, so implementing it directly
    // onto node rather than the GraphComponent trait
    pub fn update_intersection_matrix(&self, intersection_matrix: &mut IntersectionMatrix) {
        assert!(self.label.geometry_count() >= 2, "found partial label");
        intersection_matrix.set_at_least_if_valid(
            self.label.on_location(0),
            self.label.on_location(1),
            Dimensions::ZeroDimensional,
        );
        debug!(
            "updated intersection_matrix: {:?} from node: {:?}",
            intersection_matrix, self
        );
    }

    // from JTS#RelateNode, which we've squashed into the base Node class to avoid wrestling OO hierarchies into rust.
}
