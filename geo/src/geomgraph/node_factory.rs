use super::{Coordinate, Node};

// JTS: public class NodeFactory {
// JTS: /**
// JTS:  * The basic node constructor does not allow for incident edges
// JTS:  */
// JTS:   public Node createNode(Coordinate coord)
// JTS:   {
// JTS:     return new Node(coord, null);
// JTS:   }
// JTS: }

pub trait NodeFactory<F: num_traits::Float> {
    fn create_node(&self, coordinate: Coordinate<F>) -> Node<F>;
}

pub struct BasicNodeFactory;

/// The basic node constructor does not allow for incident edges
impl<F: num_traits::Float> NodeFactory<F> for BasicNodeFactory {
    fn create_node(&self, coordinate: Coordinate<F>) -> Node<F> {
        Node::new(coordinate, None)
    }
}
