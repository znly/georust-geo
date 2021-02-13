use crate::algorithm::coordinate_position::CoordinatePosition;
use crate::geomgraph::{EdgeEnd, GeometryGraph, Location, MaybeLabeledEdgeEndBundle, Position};
use crate::{Coordinate, GeoFloat, GeometryCow};
// weird circular dependency from GeomGraph to IntersectionMatrix
use crate::algorithm::dimensions::Dimensions;
use crate::algorithm::relate::IntersectionMatrix;

#[derive(Clone, Debug)]
pub(crate) struct EdgeEndBundleStar<F>
where
    F: GeoFloat,
{
    edge_map: std::collections::BTreeMap<EdgeEnd<F>, MaybeLabeledEdgeEndBundle<F>>,
    point_in_area_location: Option<[Location; 2]>,
}

impl<F> EdgeEndBundleStar<F>
where
    F: GeoFloat,
{
    pub(crate) fn new() -> Self {
        EdgeEndBundleStar {
            edge_map: std::collections::BTreeMap::new(),
            point_in_area_location: None,
        }
    }

    pub(crate) fn insert(&mut self, edge_end: EdgeEnd<F>) {
        // REVIEW: is this clone problematic? Could the coord be the `key` instead of the edge_end?
        let bundle = self
            .edge_map
            .entry(edge_end.clone())
            .or_insert(MaybeLabeledEdgeEndBundle::new(*edge_end.coordinate()));
        bundle.insert(edge_end);
    }

    pub fn update_intersection_matrix(&self, intersection_matrix: &mut IntersectionMatrix) {
        for edge_end_bundle in self.edge_end_bundles_iter() {
            edge_end_bundle.update_intersection_matrix(intersection_matrix);
            debug!(
                "updated intersection_matrix: {:?} from edge_end_bundle: {:?}",
                intersection_matrix, edge_end_bundle
            );
        }
    }
}

// From EdgeEndStar.java
// Note that JTS has EdgeEndBundleStar inherit from EdgeEndStar, but since we're only using one
// subclass (EdgeEndBundleStar) we skip the complexity of inheritance and impl this functionality
// on EdgeEndBundleStar directly. If/When we implement overlay operations we could consider
// extracting this subclass behavior
//
impl<F> EdgeEndBundleStar<F>
where
    F: GeoFloat,
{
    fn edge_end_bundles_iter(&self) -> impl Iterator<Item = &MaybeLabeledEdgeEndBundle<F>> {
        self.edge_map.values()
    }

    fn edge_end_bundles_iter_mut(
        &mut self,
    ) -> impl Iterator<Item = &mut MaybeLabeledEdgeEndBundle<F>> {
        self.edge_map.values_mut()
    }

    pub(crate) fn compute_labeling(
        &mut self,
        graph_a: &GeometryGraph<F>,
        graph_b: &GeometryGraph<F>,
    ) {
        debug!("edge_end_bundle_star: {:?}", self);
        self.compute_edge_end_labels();
        self.propagate_side_labels(0);
        self.propagate_side_labels(1);
        let mut has_dimensional_collapse_edge = [false, false];
        for edge_end in self.edge_end_bundles_iter() {
            // REVIEW: unwrap
            let label = edge_end.label().unwrap();
            for geom_index in 0..2 {
                if label.is_line(geom_index)
                    && label.on_location(geom_index) == Some(Location::Boundary)
                {
                    has_dimensional_collapse_edge[geom_index] = true;
                }
            }
        }

        // CLEANUP: in JTS this is built lazily and cached. It'll require a little borrow juggling. Is
        // that worthwhile?
        let point_in_area_location = self.point_in_area_location.or_else(|| {
            let coord = self
                .edge_end_bundles_iter()
                .next()
                .map(|edge_end_bundle| edge_end_bundle.coordinate());

            coord.map(|coord| self.build_point_in_area_location(&coord, graph_a, graph_b))
        });
        if self.point_in_area_location.is_some() && point_in_area_location.is_some() {
            self.point_in_area_location = point_in_area_location;
        }

        for edge_end_bundle in self.edge_end_bundles_iter_mut() {
            // CLEANUP: unwrap
            let label = edge_end_bundle.label_mut().unwrap();
            for geom_index in 0..2 {
                if label.is_any_empty(geom_index) {
                    let location: Location = if has_dimensional_collapse_edge[geom_index] {
                        Location::Exterior
                    } else {
                        // REVIEW: Borrowing rules forbids this.
                        // let coord = edge_end.coordinate();
                        // self.get_location(geom_index, coord, graph_a, graph_b)
                        point_in_area_location.unwrap()[geom_index]
                    };
                    label.set_all_locations_if_empty(geom_index, location);
                }
            }
        }

        debug!("edge_end_bundle_star: {:?}", self);
    }

    // TODO: support pluggable boundary node rules?
    fn compute_edge_end_labels(&mut self) {
        // REVIEW: note EdgeEndBundle is a subclass of EdgeEnd, so it's var name in JTS is a little
        // confusing
        for edge_end_bundle in self.edge_end_bundles_iter_mut() {
            edge_end_bundle.compute_label()
        }
    }

    fn build_point_in_area_location(
        &self,
        coord: &Coordinate<F>,
        graph_a: &GeometryGraph<F>,
        graph_b: &GeometryGraph<F>,
    ) -> [Location; 2] {
        [
            self.point_in_area(coord, graph_a.geometry()),
            self.point_in_area(coord, graph_b.geometry()),
        ]
    }

    fn point_in_area(&self, coord: &Coordinate<F>, geometry: &GeometryCow<F>) -> Location {
        use crate::algorithm::dimensions::HasDimensions;
        if geometry.dimensions() == Dimensions::TwoDimensional {
            geometry.coordinate_position(coord).into()
        } else {
            // if geometry is *not* an area, Coord is always Exterior
            Location::Exterior
        }
    }

    fn propagate_side_labels(&mut self, geom_index: usize) {
        let mut start_location = None;

        for edge_ends in self.edge_end_bundles_iter() {
            // REVIEW: unwrap
            let label = edge_ends.label().unwrap();
            if label.is_geom_area(geom_index) {
                if let Some(location) = label.location(geom_index, Position::Left) {
                    start_location = Some(location);
                }
            }
        }
        if start_location.is_none() {
            return;
        }
        let mut current_location = start_location.unwrap();

        for edge_ends in self.edge_end_bundles_iter_mut() {
            let label = edge_ends.label_mut().unwrap();
            if label.location(geom_index, Position::On).is_none() {
                label.set_location(geom_index, Position::On, current_location);
            }
            if label.is_geom_area(geom_index) {
                let left_location = label.location(geom_index, Position::Left);
                let right_location = label.location(geom_index, Position::Right);

                if let Some(right_location) = right_location {
                    debug_assert!(right_location == current_location, "side_location conflict with coordinate: {:?}, right_location: {:?}, current_location: {:?}", edge_ends.coordinate(), right_location, current_location);
                    assert!(left_location.is_some(), "found single null side");
                    current_location = left_location.unwrap();
                } else {
                    debug_assert!(label.location(geom_index, Position::Left).is_none());
                    label.set_location(geom_index, Position::Right, current_location);
                    label.set_location(geom_index, Position::Left, current_location);
                }
            }
        }
    }
}
