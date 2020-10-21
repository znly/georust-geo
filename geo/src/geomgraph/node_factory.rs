use super::{BasicNode, Coordinate, Node};

// JTS: public class NodeFactory {
// JTS: /**
// JTS:  * The basic node constructor does not allow for incident edges
// JTS:  */
// JTS:   public Node createNode(Coordinate coord)
// JTS:   {
// JTS:     return new Node(coord, null);
// JTS:   }
// JTS: }

pub trait NodeFactory<F: num_traits::Float, N: Node<F>> {
    fn create_node(coordinate: Coordinate<F>) -> N;
}

pub struct BasicNodeFactory;

/// The basic node constructor does not allow for incident edges
impl<F: num_traits::Float> NodeFactory<F, BasicNode<F>> for BasicNodeFactory {
    fn create_node(coordinate: Coordinate<F>) -> BasicNode<F> {
        BasicNode::new(coordinate, None)
    }
}
