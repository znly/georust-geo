use super::{Coordinate, Float, Node};

// JTS: public class NodeFactory {
// JTS: /**
// JTS:  * The basic node constructor does not allow for incident edges
// JTS:  */
// JTS:   public Node createNode(Coordinate coord)
// JTS:   {
// JTS:     return new Node(coord, null);
// JTS:   }
// JTS: }

pub(crate) trait NodeFactory<F>
where
    F: Float,
{
    fn create_node(coordinate: Coordinate<F>) -> Node<F>;
}

pub(crate) struct BasicNodeFactory;

/// The basic node constructor does not allow for incident edges
impl<F> NodeFactory<F> for BasicNodeFactory
where
    F: Float,
{
    fn create_node(coordinate: Coordinate<F>) -> Node<F> {
        Node::new(coordinate, None)
    }
}
