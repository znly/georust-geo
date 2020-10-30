use super::{BasicNode, Coordinate, Float, Node};

// JTS: public class NodeFactory {
// JTS: /**
// JTS:  * The basic node constructor does not allow for incident edges
// JTS:  */
// JTS:   public Node createNode(Coordinate coord)
// JTS:   {
// JTS:     return new Node(coord, null);
// JTS:   }
// JTS: }

pub(crate) trait NodeFactory<F, N>
where
    F: Float,
    N: Node<F>,
{
    fn create_node(coordinate: Coordinate<F>) -> N;
}

pub(crate) struct BasicNodeFactory;

/// The basic node constructor does not allow for incident edges
impl<F> NodeFactory<F, BasicNode<F>> for BasicNodeFactory
where
    F: Float,
{
    fn create_node(coordinate: Coordinate<F>) -> BasicNode<F> {
        BasicNode::new(coordinate, None)
    }
}
