use crate::geomgraph::{EdgeEndBundleStar, Float, Node, NodeFactory};
use crate::Coordinate;

// JTS: import org.locationtech.jts.geom.Coordinate;
// JTS: import org.locationtech.jts.geom.IntersectionMatrix;
// JTS: import org.locationtech.jts.geomgraph.EdgeEndStar;
// JTS: import org.locationtech.jts.geomgraph.Node;
// JTS:
// JTS: /**
// JTS:  * Represents a node in the topological graph used to compute spatial relationships.
// JTS:  *
// JTS:  * @version 1.7
// JTS:  */
// JTS: public class RelateNode
// JTS:   extends Node
// JTS: {
// JTS:
// JTS:   public RelateNode(Coordinate coord, EdgeEndStar edges)
// JTS:   {
// JTS:     super(coord, edges);
// JTS:   }
// JTS:
// JTS:   /**
// JTS:    * Update the IM with the contribution for this component.
// JTS:    * A component only contributes if it has a labelling for both parent geometries
// JTS:    */
// JTS:   protected void computeIM(IntersectionMatrix im)
// JTS:   {
// JTS:     im.setAtLeastIfValid(label.getLocation(0), label.getLocation(1), 0);
// JTS:   }
// JTS:   /**
// JTS:    * Update the IM with the contribution for the EdgeEnds incident on this node.
// JTS:    */
// JTS:   void updateIMFromEdges(IntersectionMatrix im)
// JTS:   {
// JTS:     ((EdgeEndBundleStar) edges).updateIM(im);
// JTS:   }
// JTS:
// JTS: }

// Initially I tried to represent the inheritance as a base Node trait, but it got real gnarly
// as generics permeated lots of code. Instead trying a delegate approach...
// pub(crate) struct RelateNode<F>
// where
//     F: Float,
// {
//     coordinate: Coordinate<F>,
//     label: Label,
//     edge_end_bundle_star: EdgeEndBundleStar<F>,
// }
//
// impl<F> RelateNode<F>
// where
//     F: Float,
// {
//     fn new(coordinate: Coordinate<F>, edge_end_bundle_star: EdgeEndBundleStar<F>) -> Self {
//         RelateNode {
//             coordinate,
//             edge_end_bundle_star,
//             label: Label::new_with_on_location(0, None),
//         }
//     }
//
//     pub fn set_label_on_location(&mut self, geom_index: usize, location: Location) {
//         self.label.set_on_location(geom_index, location)
//     }
//
//     pub fn edges(&self) -> &EdgeEndBundleStar<F> {
//         &self.edge_end_bundle_star
//     }
// }

pub(crate) struct RelateNodeFactory;
impl<F> NodeFactory<F> for RelateNodeFactory
where
    F: Float,
{
    fn create_node(coordinate: Coordinate<F>) -> Node<F> {
        // TODO: Add RelateNodeDelegate
        Node::new(coordinate, Some(EdgeEndBundleStar::new()))
    }
}
