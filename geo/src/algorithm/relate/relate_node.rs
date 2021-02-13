use crate::geomgraph::{EdgeEndBundleStar, Node, NodeFactory};
use crate::{Coordinate, GeoFloat};

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
    F: GeoFloat,
{
    type Edges = EdgeEndBundleStar<F>;
    fn create_node(coordinate: Coordinate<F>) -> (Node<F>, Self::Edges) {
        (Node::new(coordinate), EdgeEndBundleStar::new())
    }
}
