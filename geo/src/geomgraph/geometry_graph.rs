use super::{
    algorithm::{
        boundary_node_rule::{BoundaryNodeRule, Mod2BoundaryNodeRule},
        LineIntersector,
    },
    index::{EdgeSetIntersector, SegmentIntersector, SimpleEdgeSetIntersector},
    BasicNodeFactory, Edge, Label, Location, Node, NodeFactory, PlanarGraph, Position,
    TopologyLocation,
};

use crate::algorithm::dimensions::HasDimensions;
use crate::{Coordinate, GeoFloat, GeometryCow, Line, LineString, Point, Polygon};

use std::cell::RefCell;
use std::rc::Rc;

/// A GeometryGraph is a graph that models a given Geometry
pub(crate) struct GeometryGraph<'a, F>
where
    F: GeoFloat,
{
    arg_index: usize,
    parent_geometry: &'a GeometryCow<'a, F>,
    use_boundary_determination_rule: bool,
    planar_graph: PlanarGraph<F>,
}

// PlanarGraph delegations
// In JTS this is achieved through inheritance - GeometryGraph inherits from PlanarGraph
impl<F> GeometryGraph<'_, F>
where
    F: GeoFloat,
{
    pub fn edges(&self) -> &[Rc<RefCell<Edge<F>>>] {
        self.planar_graph.edges()
    }

    pub fn insert_edge(&mut self, edge: Edge<F>) {
        self.planar_graph.insert_edge(edge)
    }

    pub fn is_boundary_node(&self, geom_index: usize, coord: Coordinate<F>) -> bool {
        self.planar_graph.is_boundary_node(geom_index, coord)
    }

    pub fn add_node_with_coordinate(&mut self, coord: Coordinate<F>) -> &mut Node<F> {
        self.planar_graph.add_node_with_coordinate(coord)
    }

    pub fn nodes_iter(&self) -> impl Iterator<Item = &Node<F>> {
        self.planar_graph.nodes.iter().map(|t| &t.0)
    }
}

impl<'a, F> GeometryGraph<'a, F>
where
    F: GeoFloat,
{
    pub fn determine_boundary(
        boundary_node_rule: &impl BoundaryNodeRule,
        boundary_count: usize,
    ) -> Location {
        if boundary_node_rule.is_in_boundary(boundary_count) {
            Location::Boundary
        } else {
            Location::Interior
        }
    }

    // CLEANUP: hardcode the type? is perf significant?
    fn create_edge_set_intersector() -> Box<dyn EdgeSetIntersector<F>> {
        // TODO: implement more optimized edge set intersectors
        return Box::new(SimpleEdgeSetIntersector::new());
    }

    pub fn new(arg_index: usize, parent_geometry: &'a GeometryCow<F>) -> Self {
        let mut graph = GeometryGraph {
            arg_index,
            parent_geometry,
            use_boundary_determination_rule: true,
            planar_graph: PlanarGraph::new(),
        };
        graph.add_geometry(parent_geometry);
        graph
    }

    pub fn geometry(&self) -> &GeometryCow<F> {
        self.parent_geometry
    }

    fn boundary_nodes(&self) -> impl Iterator<Item = &Node<F>> {
        // TODO: should we need to memoize this like JTS does?
        self.planar_graph
            .nodes
            .boundary_nodes(self.arg_index)
            .map(|(n, _)| n)
    }

    pub fn add_geometry(&mut self, geometry: &GeometryCow<F>) {
        if geometry.is_empty() {
            return;
        }
        match geometry {
            GeometryCow::Line(line) => self.add_line(line),
            GeometryCow::Rect(rect) => {
                // PERF: avoid this conversion/clone?
                self.add_polygon(&Polygon::from(rect.clone().into_owned()).into());
            }
            GeometryCow::Point(point) => {
                self.add_point(point);
            }
            GeometryCow::Polygon(polygon) => self.add_polygon(polygon),
            GeometryCow::Triangle(triangle) => {
                // PERF: avoid this conversion/clone?
                self.add_polygon(&Polygon::from(triangle.clone().into_owned()));
            }
            GeometryCow::LineString(line_string) => self.add_line_string(line_string),
            GeometryCow::MultiPoint(multi_point) => {
                for point in &multi_point.0 {
                    self.add_point(point);
                }
            }
            GeometryCow::MultiPolygon(multi_polygon) => {
                self.use_boundary_determination_rule = false;
                for polygon in &multi_polygon.0 {
                    self.add_polygon(polygon);
                }
            }
            GeometryCow::MultiLineString(multi_line_string) => {
                for line_string in &multi_line_string.0 {
                    // PERF: can we get rid of these clones?
                    self.add_line_string(line_string);
                }
            }
            GeometryCow::GeometryCollection(geometry_collection) => {
                for geometry in geometry_collection.iter() {
                    self.add_geometry(&GeometryCow::from(geometry));
                }
            }
        }
    }

    fn add_polygon_ring(
        &mut self,
        linear_ring: &LineString<F>,
        cw_left: Location,
        cw_right: Location,
    ) {
        debug_assert!(linear_ring.is_closed());
        if linear_ring.is_empty() {
            return;
        }

        let mut coords: Vec<Coordinate<F>> = vec![];
        for coord in &linear_ring.0 {
            if coords.last() != Some(coord) {
                coords.push(*coord)
            }
        }

        if coords.len() < 4 {
            todo!("handle invalid ring")
        }
        let first_point = coords[0].clone();

        use crate::algorithm::winding_order::{Winding, WindingOrder};
        let (left, right) = if linear_ring.winding_order() == Some(WindingOrder::CounterClockwise) {
            (cw_right, cw_left)
        } else {
            (cw_left, cw_right)
        };

        let edge = Edge::new(
            coords,
            Label::new(
                self.arg_index,
                TopologyLocation::area(Location::Boundary, left, right),
            ),
        );
        // REVIEW: note we don't implement lineEdgeMap. I don't think we *need* it for the Relate
        // operations and I think it'll require moving edges into a RefCell.

        self.insert_edge(edge);

        // insert the endpoint as a node, to mark that it is on the boundary
        self.insert_point(self.arg_index, first_point, Location::Boundary);
    }

    fn add_polygon(&mut self, polygon: &Polygon<F>) {
        self.add_polygon_ring(polygon.exterior(), Location::Exterior, Location::Interior);
        for hole in polygon.interiors() {
            self.add_polygon_ring(hole, Location::Interior, Location::Exterior)
        }
    }

    fn add_line_string(&mut self, line_string: &LineString<F>) {
        let mut coords: Vec<Coordinate<F>> = vec![];
        for coord in &line_string.0 {
            if coords.last() != Some(coord) {
                coords.push(*coord)
            }
        }

        if coords.len() < 2 {
            todo!("handle invalid line string");
        }
        self.insert_boundary_point(self.arg_index, *coords.first().unwrap());
        self.insert_boundary_point(self.arg_index, *coords.last().unwrap());

        let edge = Edge::new(
            coords,
            Label::new(self.arg_index, TopologyLocation::line(Location::Interior)),
        );

        // REVIEW: note we don't implement lineEdgeMap. I don't think we *need* it for the Relate
        // operations and I think it'll require moving edges into a RefCell.

        self.insert_edge(edge);

        // REVIEW: re-ordered code to insert boundary points *before* `coords` moves into Edge::new
    }

    fn add_line(&mut self, line: &Line<F>) {
        self.insert_boundary_point(self.arg_index, line.start);
        self.insert_boundary_point(self.arg_index, line.end);

        let edge = Edge::new(
            vec![line.start, line.end],
            Label::new(self.arg_index, TopologyLocation::line(Location::Interior)),
        );

        self.insert_edge(edge);
    }

    /// Add a point computed externally.  The point is assumed to be a
    /// Point Geometry part, which has a location of INTERIOR.
    fn add_point(&mut self, point: &Point<F>) {
        self.insert_point(self.arg_index, point.clone().into(), Location::Interior);
    }

    /// Compute self-nodes, taking advantage of the Geometry type to minimize the number of
    /// intersection tests.  (E.g. rings are not tested for self-intersection, since they are
    /// assumed to be valid).
    ///
    /// @param li the LineIntersector to use
    /// @param computeRingSelfNodes if <code>false</code>, intersection checks are optimized to not test rings for self-intersection
    /// @param isDoneIfProperInt short-circuit the intersection computation if a proper intersection is found
    /// @return the computed SegmentIntersector containing information about the intersections found
    pub fn compute_self_nodes(
        &mut self,
        line_intersector: Box<dyn LineIntersector<F>>,
        compute_ring_self_nodes: bool,
        is_done_when_proper_intersection: bool,
    ) -> SegmentIntersector<F> {
        let mut segment_intersector = SegmentIntersector::new(line_intersector, true, false);

        segment_intersector.set_is_done_when_proper_intersection(is_done_when_proper_intersection);

        let mut edge_set_intersector = Self::create_edge_set_intersector();

        // optimize intersection search for valid Polygons and LinearRings
        let is_rings = match self.geometry() {
            GeometryCow::LineString(ls) => ls.is_closed(),
            GeometryCow::MultiLineString(ls) => ls.is_closed(),
            GeometryCow::Polygon(_) | GeometryCow::MultiPolygon(_) => true,
            _ => false,
        };
        let compute_all_segments = compute_ring_self_nodes || !is_rings;

        edge_set_intersector.compute_intersections(
            self.edges(),
            &mut segment_intersector,
            compute_all_segments,
        );

        // CLEANUP: read self.arg_index as property within addSelfIntersectionNodes rather than
        //          pass as param?
        self.add_self_intersection_nodes(self.arg_index);

        segment_intersector
    }

    pub fn compute_edge_intersections(
        &self,
        other: &GeometryGraph<F>,
        line_intersector: Box<dyn LineIntersector<F>>,
        include_proper: bool,
    ) -> SegmentIntersector<F> {
        let mut segment_intersector =
            SegmentIntersector::new(line_intersector, include_proper, true);
        segment_intersector.set_boundary_nodes(
            // CLEANUP: surely there's a nicer way to: `Vec<&Node> -> Vec<Node>`
            self.boundary_nodes()
                .into_iter()
                .map(|n| n.clone())
                .collect(),
            other
                .boundary_nodes()
                .into_iter()
                .map(|n| n.clone())
                .collect(),
        );

        let mut edge_set_intersector = Self::create_edge_set_intersector();
        edge_set_intersector.compute_intersections_testing_all_segments(
            self.edges(),
            other.edges(),
            &mut segment_intersector,
        );

        segment_intersector
    }

    fn insert_point(&mut self, arg_index: usize, coord: Coordinate<F>, location: Location) {
        let node: &mut Node<F> = self.add_node_with_coordinate(coord);
        node.label_mut().set_on_location(arg_index, location);
    }

    /// Adds candidate boundary points using the current {@link BoundaryNodeRule}.
    /// This is used to add the boundary points of dim-1 geometries (Curves/MultiCurves).
    fn insert_boundary_point(&mut self, arg_index: usize, coord: Coordinate<F>) {
        let node: &mut Node<F> = self.add_node_with_coordinate(coord);

        let label: &mut Label = node.label_mut();

        // the new point to insert is on a boundary
        let mut boundary_count = 1;
        // determine the current location for the point (if any)
        let location = label.location(arg_index, Position::On);
        if let Some(Location::Boundary) = location {
            boundary_count += 1;
        }

        // determine the boundary status of the point according to the Boundary Determination Rule
        // TODO: accommodate pluggable boundary node rules
        let new_location = Self::determine_boundary(&Mod2BoundaryNodeRule, boundary_count);
        label.set_on_location(arg_index, new_location);
    }
    fn add_self_intersection_nodes(&mut self, arg_index: usize) {
        let locations_and_intersections: Vec<(Location, Vec<Coordinate<F>>)> = self
            .edges()
            .into_iter()
            .map(|cell| cell.borrow())
            .map(|edge| {
                // REVIEW: move up unwrap
                let location = edge.label().on_location(arg_index).unwrap();

                let coordinates = edge
                    .edge_intersections()
                    .into_iter()
                    .map(|edge_intersection| edge_intersection.coordinate());

                (location, coordinates.collect())
            })
            .collect();

        for (location, edge_intersection_coordinates) in locations_and_intersections {
            for coordinate in edge_intersection_coordinates {
                self.add_self_intersection_node(arg_index, coordinate, location)
            }
        }
    }

    /// Add a node for a self-intersection.
    ///
    /// If the node is a potential boundary node (e.g. came from an edge which is a boundary), then
    /// insert it as a potential boundary node.  Otherwise, just add it as a regular node.
    fn add_self_intersection_node(
        &mut self,
        arg_index: usize,
        coord: Coordinate<F>,
        location: Location,
    ) {
        // if this node is already a boundary node, don't change it
        if self.is_boundary_node(arg_index, coord) {
            return;
        }

        if location == Location::Boundary && self.use_boundary_determination_rule {
            self.insert_boundary_point(arg_index, coord)
        } else {
            self.insert_point(arg_index, coord, location)
        }
    }
}
