use super::EdgeEndBundleStar;
use crate::geomgraph::{EdgeEnd, Float, GraphComponent, Label, Location, Node, NodeFactory};
use crate::Coordinate;

pub(crate) struct RelateNode<F>
where
    F: Float,
{
    coordinate: Coordinate<F>,
    label: Label,
    edge_end_bundle_star: EdgeEndBundleStar<F>,
}

impl<F> Node<F> for RelateNode<F>
where
    F: Float,
{
    // REVIEW: duplicated from BasicNode since we don't have inheritance
    fn coordinate(&self) -> &Coordinate<F> {
        &self.coordinate
    }

    // REVIEW: duplicated from BasicNode since we don't have inheritance
    // JTS:   public void add(EdgeEnd e)
    // JTS:   {
    // JTS:     // Assert: start pt of e is equal to node point
    // JTS:     edges.insert(e);
    // JTS:     e.setNode(this);
    // JTS:   }
    fn add_edge_end(&self, edge_end: EdgeEnd<F>) {
        // REVIEW: get rid of uwrap?
        // edge_end.set_node(self);
        // self.edge_end_bundle_star.insert(edge_end);
        todo!()
    }
}

impl<F> GraphComponent for RelateNode<F>
where
    F: Float,
{
    fn label(&self) -> Option<&Label> {
        Some(&self.label)
    }

    fn label_mut(&mut self) -> Option<&mut Label> {
        Some(&mut self.label)
    }

    fn set_label(&mut self, new_value: Label) {
        self.label = new_value;
    }

    fn is_isolated(&self) -> bool {
        self.label.geometry_count() == 1
    }
}

impl<F> RelateNode<F>
where
    F: Float,
{
    fn new(coordinate: Coordinate<F>, edge_end_bundle_star: EdgeEndBundleStar<F>) -> Self {
        RelateNode {
            coordinate,
            edge_end_bundle_star,
            label: Label::new_with_on_location(0, None),
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

    pub fn edges(&self) -> &EdgeEndBundleStar<F> {
        &self.edge_end_bundle_star
    }
}

pub(crate) struct RelateNodeFactory;
impl<F> NodeFactory<F, RelateNode<F>> for RelateNodeFactory
where
    F: Float,
{
    fn create_node(coordinate: Coordinate<F>) -> RelateNode<F> {
        RelateNode::new(coordinate, EdgeEndBundleStar::new())
    }
}
