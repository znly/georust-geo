use crate::geomgraph::{GraphComponent, Label, Node, NodeFactory};
use crate::Coordinate;

pub struct RelateNode<F: num_traits::Float> {
    coordinate: Coordinate<F>,
    label: Label,
    edge_end_bundle_star: EdgeEndBundleStar,
}

impl<F: num_traits::Float> Node<F> for RelateNode<F> {
    fn coordinate(&self) -> &Coordinate<F> {
        &self.coordinate
    }
}

impl<F: num_traits::Float> GraphComponent for RelateNode<F> {
    fn label(&self) -> Option<&Label> {
        Some(&self.label)
    }

    fn label_mut(&mut self) -> Option<&mut Label> {
        Some(&mut self.label)
    }

    fn set_label(&mut self, new_value: Label) {
        self.label = new_value;
    }
}

impl<F: num_traits::Float> RelateNode<F> {
    fn new(coordinate: Coordinate<F>, edge_end_bundle_star: EdgeEndBundleStar) -> Self {
        RelateNode {
            coordinate,
            edge_end_bundle_star,
            label: Label::new(0, None),
        }
    }
}

pub struct RelateNodeFactory;
impl<F: num_traits::Float> NodeFactory<F, RelateNode<F>> for RelateNodeFactory {
    fn create_node(coordinate: Coordinate<F>) -> RelateNode<F> {
        RelateNode::new(coordinate, EdgeEndBundleStar::new())
    }
}

pub struct EdgeEndBundleStar;
impl EdgeEndBundleStar {
    fn new() -> Self {
        todo!()
    }
}
