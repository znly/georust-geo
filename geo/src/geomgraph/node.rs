// JTS: import org.locationtech.jts.geom.Coordinate;
// JTS: import org.locationtech.jts.geom.IntersectionMatrix;
// JTS: import org.locationtech.jts.geom.Location;
use super::{Coordinate, EdgeEnd, EdgeEndBundleStar, Float, GraphComponent, Label, Location};

pub(crate) trait Node<F>: GraphComponent
where
    F: num_traits::Float,
{
    fn coordinate(&self) -> &Coordinate<F>;
    fn add_edge_end(&self, edge_end: EdgeEnd<F>);
}

// JTS: /**
// JTS:  * @version 1.7
// JTS:  */
// JTS: public class Node
// JTS:   extends GraphComponent
// JTS: {
#[derive(Clone)]
pub(crate) struct Node<F>
where
    F: Float,
{
    coordinate: Coordinate<F>,
    // CLEANUP: can we get rid of this Option?
    //  - BasicNodeFactory which is used by GeometryGraph via PlanarGraph passes edges: None
    //  - RelateNodeFactory which is used  by RelateComputer passes edge: EdgeEndBundleStar
    //  - TODO: in JTS the overlay operation stuff uses DirectedEdgeStar.
    // It becomes very prickly to make it generic or boxed, so sticking with just what we need for
    // relate for now.
    edges: Option<EdgeEndBundleStar<F>>,
    label: Label,
}

impl<F> GraphComponent for Node<F>
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

impl<F> Node<F>
where
    F: Float,
{
    // JTS:   protected Coordinate coord; // only non-null if this node is precise
    // JTS:   public void add(EdgeEnd e)
    // JTS:   {
    // JTS:     // Assert: start pt of e is equal to node point
    // JTS:     edges.insert(e);
    // JTS:     e.setNode(this);
    // JTS:   }
    pub fn add_edge_end(&mut self, mut edge_end: EdgeEnd<F>) {
        edge_end.set_node(self);
        // REVIEW: get rid of uwrap?
        self.edges.as_mut().unwrap().insert(edge_end);
    }

    // JTS:   protected EdgeEndStar edges;
    // JTS:
    // JTS:   public Node(Coordinate coord, EdgeEndStar edges)
    // JTS:   {
    // JTS:     this.coord = coord;
    // JTS:     this.edges = edges;
    // JTS:     label = new Label(0, Location.NONE);
    // JTS:   }
    pub fn new(coordinate: Coordinate<F>, edges: Option<EdgeEndBundleStar<F>>) -> Node<F> {
        Node {
            coordinate,
            edges,
            label: Label::new_with_geom_on_location(0, None),
        }
    }

    // JTS:   public Coordinate getCoordinate() { return coord; }
    pub fn coordinate(&self) -> &Coordinate<F> {
        &self.coordinate
    }

    // JTS:   public EdgeEndStar getEdges() { return edges; }
    // REVIEW: boxed? option?
    pub fn edges(&self) -> Option<&EdgeEndBundleStar<F>> {
        self.edges.as_ref()
    }

    pub fn edges_mut(&mut self) -> Option<&mut EdgeEndBundleStar<F>> {
        self.edges.as_mut()
    }

    // JTS:
    // JTS:   /**
    // JTS:    * Tests whether any incident edge is flagged as
    // JTS:    * being in the result.
    // JTS:    * This test can be used to determine if the node is in the result,
    // JTS:    * since if any incident edge is in the result, the node must be in the result as well.
    // JTS:    *
    // JTS:    * @return <code>true</code> if any incident edge in the in the result
    // JTS:    */
    // JTS:   public boolean isIncidentEdgeInResult()
    // JTS:   {
    // JTS:     for (Iterator it = getEdges().getEdges().iterator(); it.hasNext(); ) {
    // JTS:       DirectedEdge de = (DirectedEdge) it.next();
    // JTS:       if (de.getEdge().isInResult())
    // JTS:         return true;
    // JTS:     }
    // JTS:     return false;
    // JTS:   }
    // JTS:
    // JTS:   public boolean isIsolated()
    // JTS:   {
    // JTS:     return (label.getGeometryCount() == 1);
    // JTS:   }
    // JTS:   /**
    // JTS:    * Basic nodes do not compute IMs
    // JTS:    */
    // JTS:   protected void computeIM(IntersectionMatrix im) {}
    // JTS:   /**
    // JTS:    * Add the edge to the list of edges at this node
    // JTS:    */
    // JTS:   public void mergeLabel(Node n)
    // JTS:   {
    // JTS:     mergeLabel(n.label);
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * To merge labels for two nodes,
    // JTS:    * the merged location for each LabelElement is computed.
    // JTS:    * The location for the corresponding node LabelElement is set to the result,
    // JTS:    * as long as the location is non-null.
    // JTS:    */
    // JTS:
    // JTS:   public void mergeLabel(Label label2)
    // JTS:   {
    // JTS:     for (int i = 0; i < 2; i++) {
    // JTS:       int loc = computeMergedLocation(label2, i);
    // JTS:       int thisLoc = label.getLocation(i);
    // JTS:       if (thisLoc == Location.NONE) label.setLocation(i, loc);
    // JTS:     }
    // JTS:   }
    // JTS:
    // JTS:   public void setLabel(int argIndex, int onLocation)
    // JTS:   {
    // JTS:     if (label == null) {
    // JTS:       label = new Label(argIndex, onLocation);
    // JTS:     }
    // JTS:     else
    // JTS:       label.setLocation(argIndex, onLocation);
    // JTS:   }
    pub fn set_label_on_location(&mut self, geom_index: usize, location: Location) {
        self.label.set_on_location(geom_index, Some(location))
    }

    // JTS:   /**
    // JTS:    * Updates the label of a node to BOUNDARY,
    // JTS:    * obeying the mod-2 boundaryDetermination rule.
    // JTS:    */
    // JTS:   public void setLabelBoundary(int argIndex)
    // JTS:   {
    // JTS:     if (label == null) return;
    // JTS:
    // JTS:     // determine the current location for the point (if any)
    // JTS:     int loc = Location.NONE;
    // JTS:     if (label != null)
    // JTS:       loc = label.getLocation(argIndex);
    // JTS:     // flip the loc
    // JTS:     int newLoc;
    // JTS:     switch (loc) {
    // JTS:     case Location.BOUNDARY: newLoc = Location.INTERIOR; break;
    // JTS:     case Location.INTERIOR: newLoc = Location.BOUNDARY; break;
    // JTS:     default: newLoc = Location.BOUNDARY;  break;
    // JTS:     }
    // JTS:     label.setLocation(argIndex, newLoc);
    // JTS:   }
    pub fn set_label_boundary(&mut self, geom_index: usize) {
        let new_location = match self.label.on_location(geom_index) {
            Some(Location::Boundary) => Location::Interior,
            Some(Location::Interior) => Location::Boundary,
            None | Some(Location::Exterior) => Location::Boundary,
        };
        self.label.set_on_location(geom_index, Some(new_location));
    }

    // JTS:
    // JTS:   /**
    // JTS:    * The location for a given eltIndex for a node will be one
    // JTS:    * of { null, INTERIOR, BOUNDARY }.
    // JTS:    * A node may be on both the boundary and the interior of a geometry;
    // JTS:    * in this case, the rule is that the node is considered to be in the boundary.
    // JTS:    * The merged location is the maximum of the two input values.
    // JTS:    */
    // JTS:   int computeMergedLocation(Label label2, int eltIndex)
    // JTS:   {
    // JTS:     int loc = Location.NONE;
    // JTS:     loc = label.getLocation(eltIndex);
    // JTS:     if (! label2.isNull(eltIndex)) {
    // JTS:         int nLoc = label2.getLocation(eltIndex);
    // JTS:         if (loc != Location.BOUNDARY) loc = nLoc;
    // JTS:     }
    // JTS:     return loc;
    // JTS:   }
    // JTS:
    // JTS:   public void print(PrintStream out)
    // JTS:   {
    // JTS:     out.println("node " + coord + " lbl: " + label);
    // JTS:   }
    // JTS: }
}
