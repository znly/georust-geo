mod edge;
mod edge_intersection;
mod edge_intersection_list;
mod geometry_graph;
mod index;
mod line_intersector;
mod robust_line_intersector;

use edge::Edge;
use edge_intersection::EdgeIntersection;
use edge_intersection_list::EdgeIntersectionList;
use geometry_graph::GeometryGraph;
use line_intersector::{Intersection, LineIntersector};
use robust_line_intersector::RobustLineIntersector;

use geo_types::Coordinate;
// TODO move to own file
#[derive(PartialEq)]
struct Node<F: num_traits::Float> {
    coordinate: Coordinate<F>,
}

impl<F: num_traits::Float> Node<F> {
    fn get_coordinate(&self) -> Coordinate<F> {
        self.coordinate
    }
}
