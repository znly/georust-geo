use super::{
    algorithm::boundary_node_rule::Mod2BoundaryNodeRule, Edge, EdgeEnd, GeometryGraph, Label,
    Location, Position,
};
use crate::{Coordinate, GeoFloat};
// weird circular dependency from GeomGraph to IntersectionMatrix
use crate::algorithm::relate::IntersectionMatrix;

// NOTE: in JTS this is in the algorithm::relate package, but because
// we moved EdgeEndBundleStar into geomgraph, we also must move EdgeEndBundle
// (we moved EdgeEndBundleStar rather than wrestle with inheritance of EdgeEndBundle)

#[derive(Clone, Debug)]
pub(crate) struct MaybeLabeledEdgeEndBundle<F>
where
    F: GeoFloat,
{
    coordinate: Coordinate<F>,
    label: Option<Label>,
    edge_ends: Vec<EdgeEnd<F>>,
}

impl<F> MaybeLabeledEdgeEndBundle<F>
where
    F: GeoFloat,
{
    pub(crate) fn new(coordinate: Coordinate<F>) -> Self {
        // REVIEW: there's a lot more to the JTS initializer since EdgeEndBundle inherits from
        // EdgeEnd, but let's see if it's necessary.
        Self {
            coordinate,
            edge_ends: vec![],
            label: None,
        }
    }

    pub fn label(&self) -> Option<&Label> {
        self.label.as_ref()
    }
    pub fn label_mut(&mut self) -> Option<&mut Label> {
        self.label.as_mut()
    }

    // CLEANUP: borrow?
    pub fn coordinate(&self) -> Coordinate<F> {
        self.coordinate
    }

    fn edge_ends_iter(&self) -> impl Iterator<Item = &EdgeEnd<F>> {
        self.edge_ends.iter()
    }

    fn edge_ends_iter_mut(&mut self) -> impl Iterator<Item = &mut EdgeEnd<F>> {
        self.edge_ends.iter_mut()
    }

    pub(crate) fn insert(&mut self, edge_end: EdgeEnd<F>) {
        self.edge_ends.push(edge_end);
    }

    // TODO: support pluggable boundary node rule?
    pub(crate) fn compute_label(&mut self) {
        let is_area = self
            .edge_ends_iter()
            .any(|edge_end| edge_end.label().is_area());

        if is_area {
            self.label = Some(Label::empty_area());
        } else {
            self.label = Some(Label::empty_line());
        }

        for i in 0..2 {
            self.compute_label_on(i);
            if is_area {
                self.compute_label_sides(i);
            }
        }
    }

    fn compute_label_on(&mut self, geom_index: usize) {
        let mut boundary_count = 0;
        let mut found_interior = false;

        for edge_end in self.edge_ends_iter() {
            match edge_end.label().on_location(geom_index) {
                Some(Location::Boundary) => {
                    boundary_count += 1;
                }
                Some(Location::Interior) => {
                    found_interior = true;
                }
                None | Some(Location::Exterior) => {}
            }
        }

        let mut location = None;
        if found_interior {
            location = Some(Location::Interior);
        }

        // TEST: Is it possible to be in interior *and* boundary? Is it OK that we "prefer"
        // boundary in that case?
        if boundary_count > 0 {
            location = Some(GeometryGraph::<'_, F>::determine_boundary(
                &Mod2BoundaryNodeRule,
                boundary_count,
            ));
        }

        if let Some(location) = location {
            self.label
                .as_mut()
                .map(|l| l.set_on_location(geom_index, location));
        } else {
            self.label.as_ref().map(|l| {
                // This is technically a diversion from JTS, but I don't think we'd ever
                // get here, unless `l.on_location` was *already* None, in which cases this is a
                // no-op, so assert that assumption.
                // If this assert is rightfully triggered, we may need to add a method like
                // `l.clear_on_location(geom_index)`
                debug_assert!(
                    l.on_location(geom_index).is_none(),
                    "diverging from JTS, which would have replaced the existing Location with None"
                );
            });
        }
    }

    fn compute_label_sides(&mut self, geom_index: usize) {
        self.compute_label_side(geom_index, Position::Left);
        self.compute_label_side(geom_index, Position::Right);
    }

    fn compute_label_side(&mut self, geom_index: usize, side: Position) {
        let mut location = None;
        for edge_end in self.edge_ends_iter_mut() {
            if edge_end.label().is_area() {
                match edge_end.label_mut().location(geom_index, side) {
                    Some(Location::Interior) => {
                        location = Some(Location::Interior);
                    }
                    Some(Location::Exterior) => {
                        location = Some(Location::Exterior);
                    }
                    None | Some(Location::Boundary) => {}
                }
            }
        }

        if let Some(location) = location {
            self.label
                .as_mut()
                .map(|l| l.set_location(geom_index, side, location));
        }
    }

    pub fn update_intersection_matrix(&self, intersection_matrix: &mut IntersectionMatrix) {
        // REVIEW: unwrap
        Edge::<F>::update_intersection_matrix(self.label().unwrap(), intersection_matrix);
    }
}
