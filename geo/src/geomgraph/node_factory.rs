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
