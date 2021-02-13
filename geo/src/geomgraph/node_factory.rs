use super::Node;
use crate::{Coordinate, GeoFloat};

pub(crate) trait NodeFactory<F>
where
    F: GeoFloat,
{
    type Edges;
    fn create_node(coordinate: Coordinate<F>) -> (Node<F>, Self::Edges);
}

pub(crate) struct BasicNodeFactory;

/// The basic node constructor does not allow for incident edges
impl<F> NodeFactory<F> for BasicNodeFactory
where
    F: GeoFloat,
{
    type Edges = ();
    fn create_node(coordinate: Coordinate<F>) -> (Node<F>, Self::Edges) {
        (Node::new(coordinate), ())
    }
}
