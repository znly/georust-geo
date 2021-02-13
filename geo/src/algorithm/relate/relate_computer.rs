use super::{EdgeEndBuilder, IntersectionMatrix, RelateNodeFactory};
use crate::algorithm::dimensions::{Dimensions, HasDimensions};
use crate::geomgraph::{
    algorithm::RobustLineIntersector, index::SegmentIntersector, Edge, EdgeEnd, GeometryGraph,
    Location, Node, NodeMap,
};
use crate::{GeoFloat, GeometryCow};

use std::cell::RefCell;
use std::rc::Rc;

pub(crate) struct RelateComputer<'a, F>
where
    F: GeoFloat,
{
    graph_a: GeometryGraph<'a, F>,
    graph_b: GeometryGraph<'a, F>,
    nodes: NodeMap<F, RelateNodeFactory>,
    line_intersector: RobustLineIntersector<F>,
    isolated_edges: Vec<Rc<RefCell<Edge<F>>>>,
}

impl<'a, F> RelateComputer<'a, F>
where
    F: 'static + GeoFloat,
{
    pub fn new(
        geom_a: &'a GeometryCow<'a, F>,
        geom_b: &'a GeometryCow<'a, F>,
    ) -> RelateComputer<'a, F> {
        Self {
            graph_a: GeometryGraph::new(0, geom_a),
            graph_b: GeometryGraph::new(1, geom_b),
            nodes: NodeMap::new(),
            isolated_edges: vec![],
            line_intersector: RobustLineIntersector::new(),
        }
    }

    pub fn compute_intersection_matrix(&mut self) -> IntersectionMatrix {
        let mut intersection_matrix = IntersectionMatrix::empty();
        // since Geometries are finite and embedded in a 2-D space, the EE element must always be 2
        intersection_matrix.set(
            Location::Exterior,
            Location::Exterior,
            Dimensions::TwoDimensional,
        );

        // if the Geometries don't overlap, we can skip most of the work
        use crate::algorithm::bounding_rect::BoundingRect;
        match (
            self.graph_a.geometry().bounding_rect(),
            self.graph_b.geometry().bounding_rect(),
        ) {
            (Some(bounding_rect_a), Some(bounding_rect_b)) => {
                use crate::algorithm::intersects::Intersects;
                if !bounding_rect_a.intersects(&bounding_rect_b) {
                    self.compute_disjoint_intersection_matrix(&mut intersection_matrix);
                    return intersection_matrix;
                }
            }
            _ => {
                self.compute_disjoint_intersection_matrix(&mut intersection_matrix);
                return intersection_matrix;
            }
        }

        // REVIEW: In JTS, second `false` is implied via a default arg from an overload
        // REVIEW: In JTS, self.line_intersection is just passed in. But it's mutated - seems like
        //         a bad idea and rust won't allow it. So we clone.
        self.graph_a
            .compute_self_nodes(Box::new(self.line_intersector.clone()), false, false);
        self.graph_b
            .compute_self_nodes(Box::new(self.line_intersector.clone()), false, false);

        // compute intersections between edges of the two input geometries
        let segment_intersector = self.graph_a.compute_edge_intersections(
            &self.graph_b,
            Box::new(self.line_intersector.clone()),
            false,
        );

        self.compute_intersection_nodes(0);
        self.compute_intersection_nodes(1);
        // Copy the labelling for the nodes in the parent Geometries.  These override any labels
        // determined by intersections between the geometries.
        self.copy_nodes_and_labels(0);
        self.copy_nodes_and_labels(1);
        // complete the labelling for any nodes which only have a label for a single geometry
        self.label_isolated_nodes();
        // If a proper intersection was found, we can set a lower bound on the IM.
        self.compute_proper_intersection_im(&segment_intersector, &mut intersection_matrix);
        // Now process improper intersections
        // (eg where one or other of the geometries has a vertex at the intersection point)
        // We need to compute the edge graph at all nodes to determine the IM.
        let edge_end_builder = EdgeEndBuilder::new();
        let edge_ends_a: Vec<_> = edge_end_builder.compute_ends_for_edges(self.graph_a.edges());
        // Fails - len() == 6, which is inconsistent with JTS
        // assert_eq!(edge_ends_a.len(), 2);
        self.insert_edge_ends(edge_ends_a);
        let edge_ends_b: Vec<_> = edge_end_builder.compute_ends_for_edges(self.graph_b.edges());
        self.insert_edge_ends(edge_ends_b);
        self.label_node_edges();

        self.label_isolated_edges(0, 1);
        self.label_isolated_edges(1, 0);

        debug!(
            "before update_intersection_matrix: {:?}",
            &intersection_matrix
        );
        self.update_intersection_matrix(&mut intersection_matrix);

        intersection_matrix
    }

    fn insert_edge_ends(&mut self, edge_ends: Vec<EdgeEnd<F>>) {
        for edge_end in edge_ends {
            let (_node, edges) = self.nodes.add_node_with_coordinate(*edge_end.coordinate());
            edges.insert(edge_end);
        }
    }

    fn compute_proper_intersection_im(
        &mut self,
        segment_intersector: &SegmentIntersector<F>,
        intersection_matrix: &mut IntersectionMatrix,
    ) {
        // If a proper intersection is found, we can set a lower bound on the IM.
        let dim_a = self.graph_a.geometry().dimensions();
        let dim_b = self.graph_b.geometry().dimensions();

        let has_proper = segment_intersector.has_proper_intersection();
        let has_proper_interior = segment_intersector.has_proper_interior_intersection();

        debug_assert!(
            (dim_a != Dimensions::ZeroDimensional && dim_b != Dimensions::ZeroDimensional)
                || (!has_proper && !has_proper_interior)
        );

        match (dim_a, dim_b) {
            // If edge segments of Areas properly intersect, the areas must properly overlap.
            (Dimensions::TwoDimensional, Dimensions::TwoDimensional) => {
                if has_proper {
                    intersection_matrix.set_at_least_from_string("212101212");
                }
            }

            // If a Line segment properly intersects an edge segment of an Area, it follows that
            // the Interior of the Line intersects the Boundary of the Area.
            // If the intersection is a proper *interior* intersection, then there is an
            // Interior-Interior intersection too.
            // Note that it does not follow that the Interior of the Line intersects the Exterior
            // of the Area, since there may be another Area component which contains the rest of the Line.
            (Dimensions::TwoDimensional, Dimensions::OneDimensional) => {
                if has_proper {
                    intersection_matrix.set_at_least_from_string("FFF0FFFF2");
                }

                if has_proper_interior {
                    intersection_matrix.set_at_least_from_string("1FFFFF1FF");
                }
            }

            (Dimensions::OneDimensional, Dimensions::TwoDimensional) => {
                if has_proper {
                    intersection_matrix.set_at_least_from_string("F0FFFFFF2");
                }

                if has_proper_interior {
                    intersection_matrix.set_at_least_from_string("1F1FFFFFF");
                }
            }

            // If edges of LineStrings properly intersect *in an interior point*, all we can deduce
            // is that the interiors intersect.  (We can NOT deduce that the exteriors intersect,
            // since some other segments in the geometries might cover the points in the
            // neighbourhood of the intersection.)
            // It is important that the point be known to be an interior point of both Geometries,
            // since it is possible in a self-intersecting geometry to have a proper intersection
            // on one segment that is also a boundary point of another segment.
            (Dimensions::OneDimensional, Dimensions::OneDimensional) => {
                if has_proper_interior {
                    intersection_matrix.set_at_least_from_string("0FFFFFFFF");
                }
            }
            _ => {}
        }
    }

    /// Copy all nodes from an arg geometry into this graph.
    ///
    /// The node label in the arg geometry overrides any previously computed label for that
    /// argIndex.  (E.g. a node may be an intersection node with a computed label of BOUNDARY, but
    /// in the original arg Geometry it is actually in the interior due to the Boundary
    /// Determination Rule)
    fn copy_nodes_and_labels(&mut self, geom_index: usize) {
        let graph = if geom_index == 0 {
            &self.graph_a
        } else {
            assert!(geom_index == 1);
            &self.graph_b
        };
        for graph_node in graph.nodes_iter() {
            let new_node = self
                .nodes
                .add_node_with_coordinate(*graph_node.coordinate());
            // CLEANUP: on_location().unwrap - can we get rid of it or check for it?
            new_node.0.set_label_on_location(
                geom_index,
                graph_node.label().on_location(geom_index).unwrap(),
            );
        }
    }

    /// Insert nodes for all intersections on the edges of a Geometry.  
    ///
    /// Label the created nodes the same as the edge label if they do not already have a label.
    /// This allows nodes created by either self-intersections or mutual intersections to be
    /// labelled.  
    ///
    /// Endpoint nodes will already be labeled from when they were inserted.
    fn compute_intersection_nodes(&mut self, geom_index: usize) {
        let graph = if geom_index == 0 {
            &self.graph_a
        } else {
            assert!(geom_index == 1);
            &self.graph_b
        };

        for edge in graph.edges() {
            let edge = edge.borrow();

            let edge_location = edge.label().on_location(geom_index);
            for edge_intersection in edge.edge_intersections() {
                let (new_node, _edges) = self
                    .nodes
                    .add_node_with_coordinate(edge_intersection.coordinate());

                if edge_location == Some(Location::Boundary) {
                    new_node.set_label_boundary(geom_index);
                } else {
                    if new_node.label().is_empty(geom_index) {
                        new_node.set_label_on_location(geom_index, Location::Interior);
                    }
                }
            }
        }
    }

    /// If the Geometries are disjoint, we need to enter their dimension and boundary dimension in
    /// the Ext rows in the IM
    fn compute_disjoint_intersection_matrix(&self, intersection_matrix: &mut IntersectionMatrix) {
        {
            let geometry_a = self.graph_a.geometry();
            let dimensions = geometry_a.dimensions();
            if dimensions != Dimensions::Empty {
                intersection_matrix.set(Location::Interior, Location::Exterior, dimensions);

                let boundary_dimensions = geometry_a.boundary_dimensions();
                if boundary_dimensions != Dimensions::Empty {
                    intersection_matrix.set(
                        Location::Boundary,
                        Location::Exterior,
                        boundary_dimensions,
                    );
                }
            }
        }

        {
            let geometry_b = self.graph_b.geometry();
            let dimensions = geometry_b.dimensions();
            if dimensions != Dimensions::Empty {
                intersection_matrix.set(Location::Exterior, Location::Interior, dimensions);

                let boundary_dimensions = geometry_b.boundary_dimensions();
                if boundary_dimensions != Dimensions::Empty {
                    intersection_matrix.set(
                        Location::Exterior,
                        Location::Boundary,
                        boundary_dimensions,
                    );
                }
            }
        }
    }

    fn label_node_edges(&mut self) {
        for (_node, edges) in self.nodes.iter_mut() {
            edges.compute_labeling(&self.graph_a, &self.graph_b);
        }
    }

    fn update_intersection_matrix(&self, intersection_matrix: &mut IntersectionMatrix) {
        debug!(
            "before updated_intersection_matrix(isolated_edges): {:?}",
            intersection_matrix
        );
        for isolated_edge in &self.isolated_edges {
            let edge = isolated_edge.borrow();
            Edge::<F>::update_intersection_matrix(edge.label(), intersection_matrix);
            debug!(
                "after updated_intersection_matrix(isolated_edge: {:?}, label: {:?}): {:?}",
                edge,
                edge.label(),
                intersection_matrix
            );
        }

        for (node, edges) in self.nodes.iter() {
            node.update_intersection_matrix(intersection_matrix);
            edges.update_intersection_matrix(intersection_matrix);
        }
    }

    fn label_isolated_edges(&mut self, this_index: usize, target_index: usize) {
        let (this_graph, target_graph) = if this_index == 0 {
            (&self.graph_a, &self.graph_b)
        } else {
            (&self.graph_b, &self.graph_a)
        };

        for edge in this_graph.edges() {
            let mut mut_edge = edge.borrow_mut();
            if mut_edge.is_isolated() {
                Self::label_isolated_edge(&mut mut_edge, target_index, target_graph.geometry());
                self.isolated_edges.push(edge.clone());
            }
        }
    }

    fn label_isolated_edge(edge: &mut Edge<F>, target_index: usize, target: &GeometryCow<F>) {
        if target.dimensions() > Dimensions::ZeroDimensional {
            // REVIEW: unwrap
            use crate::algorithm::coordinate_position::CoordinatePosition;
            let location = target
                .coordinate_position(edge.coordinate().unwrap())
                .into();

            edge.label_mut().set_all_locations(target_index, location);
        } else {
            edge.label_mut()
                .set_all_locations(target_index, Location::Exterior);
        }
    }

    /// Isolated nodes are nodes whose labels are incomplete (e.g. the location for one Geometry is
    /// null).  
    /// This is the case because nodes in one graph which don't intersect nodes in the other
    /// are not completely labelled by the initial process of adding nodes to the nodeList.  To
    /// complete the labelling we need to check for nodes that lie in the interior of edges, and in
    /// the interior of areas.
    fn label_isolated_nodes(&mut self) {
        let geometry_a = self.graph_a.geometry();
        let geometry_b = self.graph_b.geometry();
        for (node, _edges) in self.nodes.iter_mut() {
            let label = node.label();
            // isolated nodes should always have at least one geometry in their label
            debug_assert!(label.geometry_count() > 0, "node with empty label found");
            if node.is_isolated() {
                if label.is_empty(0) {
                    Self::label_isolated_node(node, 0, geometry_a)
                } else {
                    Self::label_isolated_node(node, 1, geometry_b)
                }
            }
        }
    }

    fn label_isolated_node(node: &mut Node<F>, target_index: usize, geometry: &GeometryCow<F>) {
        use crate::algorithm::coordinate_position::CoordinatePosition;
        let location = geometry.coordinate_position(node.coordinate()).into();
        node.label_mut().set_all_locations(target_index, location);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use geo_types::{polygon, Geometry};

    #[test]
    fn test_disjoint() {
        let square_a: Geometry<f64> = polygon![
            (x: 0., y: 0.),
            (x: 0., y: 20.),
            (x: 20., y: 20.),
            (x: 20., y: 0.),
            (x: 0., y: 0.),
        ]
        .into();

        let square_b: Geometry<f64> = polygon![
            (x: 55., y: 55.),
            (x: 50., y: 60.),
            (x: 60., y: 60.),
            (x: 60., y: 55.),
            (x: 55., y: 55.),
        ]
        .into();

        let gc1 = GeometryCow::from(&square_a);
        let gc2 = GeometryCow::from(&square_b);
        let mut relate_computer = RelateComputer::new(&gc1, &gc2);
        let intersection_matrix = relate_computer.compute_intersection_matrix();
        assert_eq!(
            intersection_matrix,
            IntersectionMatrix::from_str("FF2FF1212")
        );
    }

    #[test]
    fn test_a_contains_b() {
        let square_a: Geometry<f64> = polygon![
            (x: 0., y: 0.),
            (x: 0., y: 20.),
            (x: 20., y: 20.),
            (x: 20., y: 0.),
            (x: 0., y: 0.),
        ]
        .into();

        let square_b: Geometry<f64> = polygon![
            (x: 5., y: 5.),
            (x: 5., y: 10.),
            (x: 10., y: 10.),
            (x: 10., y: 5.),
            (x: 5., y: 5.),
        ]
        .into();

        let gca = GeometryCow::from(&square_a);
        let gcb = GeometryCow::from(&square_b);
        let mut relate_computer = RelateComputer::new(&gca, &gcb);
        let intersection_matrix = relate_computer.compute_intersection_matrix();
        assert_eq!(
            intersection_matrix,
            IntersectionMatrix::from_str("212FF1FF2")
        );
    }

    #[test]
    fn test_a_overlaps_b() {
        let square_a: Geometry<f64> = polygon![
            (x: 0., y: 0.),
            (x: 0., y: 20.),
            (x: 20., y: 20.),
            (x: 20., y: 0.),
            (x: 0., y: 0.),
        ]
        .into();

        let square_b: Geometry<f64> = polygon![
            (x: 5., y: 5.),
            (x: 5., y: 30.),
            (x: 30., y: 30.),
            (x: 30., y: 5.),
            (x: 5., y: 5.),
        ]
        .into();

        let gca = &GeometryCow::from(&square_a);
        let gcb = &GeometryCow::from(&square_b);
        let mut relate_computer = RelateComputer::new(gca, gcb);
        let intersection_matrix = relate_computer.compute_intersection_matrix();
        assert_eq!(
            intersection_matrix,
            IntersectionMatrix::from_str("212101212")
        );
    }
}
