#![allow(dead_code)]
#![allow(unused_imports)]

mod edge;
mod edge_intersection;
mod edge_intersection_list;
mod geometry_graph;
mod graph_component;
mod index;
mod label;
mod node;
mod node_factory;
mod node_map;
mod planar_graph;
mod topology_location;

use edge::Edge;
use edge_intersection::EdgeIntersection;
use edge_intersection_list::EdgeIntersectionList;
pub use geometry_graph::GeometryGraph;
pub use graph_component::GraphComponent;
pub(crate) use label::Label;
pub(crate) use node::{BasicNode, Node};
pub(crate) use node_factory::{BasicNodeFactory, NodeFactory};
pub use node_map::NodeMap;
use planar_graph::PlanarGraph;
use topology_location::TopologyLocation;

// in just algorithm is an external pacakage (top level, still part of JTS) - not witin geomgraph
pub(crate) mod algorithm;

use geo_types::{Coordinate, Geometry};

// CLEANUP: use geo::kernels::Orientation instead?
#[derive(Copy, Clone, PartialEq, Eq)]
pub enum Position {
    // CLEANUP: get rid of the explicit discrimanator?
    On = 0,
    Left = 1,
    Right = 2,
}

// CLEANUP: geo::utils::CoordPos instead?
#[derive(Copy, Clone, PartialEq, Eq)]
pub enum Location {
    // CLEANUP: get rid of the explicit discrimanator?
    Interior = 0,
    Boundary = 1,
    Exterior = 2,
}

/// TODO: actually implement EdgeEndStar for incident edge intersections
#[derive(Clone)]
pub struct EdgeEndStar;
