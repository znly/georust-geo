use super::{BasicNodeFactory, Edge, Label, Location, Node, NodeMap};
use crate::{Coordinate, GeoFloat};

use std::cell::RefCell;
use std::rc::Rc;

/// The computation of the [IntersectionMatrix] relies on the use of a
/// structure called a "topology graph". The topology graph contains nodes and
/// edges corresponding to the nodes and line segments of a [Geometry]. Each
/// node and edge in the graph is labeled with its topological location
/// relative to the source geometry.
///
/// Note that there is no requirement that points of self-intersection be a
/// vertex.  Thus to obtain a correct topology graph, [Geometry] must be
/// self-noded before constructing their graphs.
///
/// Two fundamental operations are supported by topology graphs:
///   - Computing the intersections between all the edges and nodes of a single graph
///   - Computing the intersections between the edges and nodes of two different graphs
pub(crate) struct PlanarGraph<F: GeoFloat> {
    pub(crate) nodes: NodeMap<F, BasicNodeFactory>,
    edges: Vec<Rc<RefCell<Edge<F>>>>,
}

impl<F: GeoFloat> PlanarGraph<F> {
    pub fn edges(&self) -> &[Rc<RefCell<Edge<F>>>] {
        &self.edges
    }

    pub fn new() -> Self {
        PlanarGraph {
            nodes: NodeMap::new(),
            edges: vec![],
        }
    }

    pub fn is_boundary_node(&self, geom_index: usize, coord: Coordinate<F>) -> bool {
        self.nodes
            .find(coord)
            .and_then(|(node, _)| node.label().on_location(geom_index))
            .map(|location| location == Location::Boundary)
            .unwrap_or(false)
    }

    pub fn insert_edge(&mut self, edge: Edge<F>) {
        self.edges.push(Rc::new(RefCell::new(edge)));
    }

    pub fn add_node_with_coordinate(&mut self, coord: Coordinate<F>) -> &mut Node<F> {
        &mut self.nodes.add_node_with_coordinate(coord).0
    }
}
