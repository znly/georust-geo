use crate::geomgraph::{GraphComponent, Label, Location, Node, NodeFactory};
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

    pub fn set_label_boundary(&mut self, geom_index: usize) {
        // REVIEW: This is actually from JTS's base Node class, but so far only used in RelateNode
        // JTS: /**
        // JTS:  * Updates the label of a node to BOUNDARY,
        // JTS:  * obeying the mod-2 boundaryDetermination rule.
        // JTS:  */
        // JTS: public void setLabelBoundary(int argIndex)
        // JTS: {
        // JTS:   if (label == null) return;

        // JTS:   // determine the current location for the point (if any)
        // JTS:   int loc = Location.NONE;
        // JTS:   if (label != null)
        // JTS:     loc = label.getLocation(argIndex);
        // JTS:   // flip the loc
        // JTS:   int newLoc;
        // JTS:   switch (loc) {
        // JTS:   case Location.BOUNDARY: newLoc = Location.INTERIOR; break;
        // JTS:   case Location.INTERIOR: newLoc = Location.BOUNDARY; break;
        // JTS:   default: newLoc = Location.BOUNDARY;  break;
        // JTS:   }
        // JTS:   label.setLocation(argIndex, newLoc);
        // JTS: }
        // this might be an overzealous assert - JTS doesn't leverage Optional types
        debug_assert!(self.label.on_location(geom_index).is_some());
        let new_location = match self.label.on_location(geom_index) {
            Some(Location::Boundary) => Location::Interior,
            Some(Location::Interior) => Location::Boundary,
            None | Some(Location::Exterior) => Location::Boundary,
        };
        self.label.set_on_location(geom_index, new_location);
    }

    pub fn set_label_on_location(&mut self, geom_index: usize, location: Location) {
        self.label.set_on_location(geom_index, location)
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
