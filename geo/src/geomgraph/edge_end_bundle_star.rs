use crate::geomgraph::{EdgeEnd, EdgeEndBundle, Float, GeometryGraph};

// JTS: /**
// JTS:  * An ordered list of {@link EdgeEndBundle}s around a {@link RelateNode}.
// JTS:  * They are maintained in CCW order (starting with the positive x-axis) around the node
// JTS:  * for efficient lookup and topology building.
// JTS:  *
// JTS:  * @version 1.7
// JTS:  */
// JTS: public class EdgeEndBundleStar
// JTS:   extends EdgeEndStar
// JTS: {
// JTS:   /**
// JTS:    * Creates a new empty EdgeEndBundleStar
// JTS:    */
// JTS:   public EdgeEndBundleStar() {
// JTS:   }
// JTS:
#[derive(Clone)]
pub(crate) struct EdgeEndBundleStar<F>
where
    F: Float,
{
    edge_map: std::collections::BTreeMap<EdgeEnd<F>, EdgeEndBundle<F>>,
}

impl<F> EdgeEndBundleStar<F>
where
    F: Float,
{
    pub(crate) fn new() -> Self {
        EdgeEndBundleStar {
            edge_map: std::collections::BTreeMap::new(),
        }
    }

    // JTS:   /**
    // JTS:    * Insert a EdgeEnd in order in the list.
    // JTS:    * If there is an existing EdgeStubBundle which is parallel, the EdgeEnd is
    // JTS:    * added to the bundle.  Otherwise, a new EdgeEndBundle is created
    // JTS:    * to contain the EdgeEnd.
    // JTS:    * <br>
    // JTS:    */
    // JTS:   public void insert(EdgeEnd e)
    // JTS:   {
    // JTS:     EdgeEndBundle eb = (EdgeEndBundle) edgeMap.get(e);
    // JTS:     if (eb == null) {
    // JTS:       eb = new EdgeEndBundle(e);
    // JTS:       insertEdgeEnd(e, eb);
    // JTS:     }
    // JTS:     else {
    // JTS:       eb.insert(e);
    // JTS:     }
    // JTS:   }
    pub(crate) fn insert(&mut self, edge_end: EdgeEnd<F>) {
        let bundle = self
            .edge_map
            .entry(edge_end.clone())
            .or_insert(EdgeEndBundle::new());
        bundle.insert(edge_end);
    }

    // JTS:   /**
    // JTS:    * Update the IM with the contribution for the EdgeStubs around the node.
    // JTS:    */
    // JTS:   void updateIM(IntersectionMatrix im)
    // JTS:   {
    // JTS:     for (Iterator it = iterator(); it.hasNext(); ) {
    // JTS:       EdgeEndBundle esb = (EdgeEndBundle) it.next();
    // JTS:       esb.updateIM(im);
    // JTS:     }
    // JTS:   }
    // JTS:
    // JTS: }
    //
    //
    // From EdgeEndStar base class in JTS, which we don't implement
    // JTS:   public void computeLabelling(GeometryGraph[] geomGraph)
    // JTS:   {
    // JTS:     computeEdgeEndLabels(geomGraph[0].getBoundaryNodeRule());
    // JTS:     // Propagate side labels  around the edges in the star
    // JTS:     // for each parent Geometry
    // JTS: //Debug.print(this);
    // JTS:     propagateSideLabels(0);
    // JTS: //Debug.print(this);
    // JTS: //Debug.printIfWatch(this);
    // JTS:     propagateSideLabels(1);
    // JTS: //Debug.print(this);
    // JTS: //Debug.printIfWatch(this);
    // JTS:
    // JTS:     /**
    // JTS:      * If there are edges that still have null labels for a geometry
    // JTS:      * this must be because there are no area edges for that geometry incident on this node.
    // JTS:      * In this case, to label the edge for that geometry we must test whether the
    // JTS:      * edge is in the interior of the geometry.
    // JTS:      * To do this it suffices to determine whether the node for the edge is in the interior of an area.
    // JTS:      * If so, the edge has location INTERIOR for the geometry.
    // JTS:      * In all other cases (e.g. the node is on a line, on a point, or not on the geometry at all) the edge
    // JTS:      * has the location EXTERIOR for the geometry.
    // JTS:      * <p>
    // JTS:      * Note that the edge cannot be on the BOUNDARY of the geometry, since then
    // JTS:      * there would have been a parallel edge from the Geometry at this node also labelled BOUNDARY
    // JTS:      * and this edge would have been labelled in the previous step.
    // JTS:      * <p>
    // JTS:      * This code causes a problem when dimensional collapses are present, since it may try and
    // JTS:      * determine the location of a node where a dimensional collapse has occurred.
    // JTS:      * The point should be considered to be on the EXTERIOR
    // JTS:      * of the polygon, but locate() will return INTERIOR, since it is passed
    // JTS:      * the original Geometry, not the collapsed version.
    // JTS:      *
    // JTS:      * If there are incident edges which are Line edges labelled BOUNDARY,
    // JTS:      * then they must be edges resulting from dimensional collapses.
    // JTS:      * In this case the other edges can be labelled EXTERIOR for this Geometry.
    // JTS:      *
    // JTS:      * MD 8/11/01 - NOT TRUE!  The collapsed edges may in fact be in the interior of the Geometry,
    // JTS:      * which means the other edges should be labelled INTERIOR for this Geometry.
    // JTS:      * Not sure how solve this...  Possibly labelling needs to be split into several phases:
    // JTS:      * area label propagation, symLabel merging, then finally null label resolution.
    // JTS:      */
    // JTS:     boolean[] hasDimensionalCollapseEdge = { false, false };
    // JTS:     for (Iterator it = iterator(); it.hasNext(); ) {
    // JTS:       EdgeEnd e = (EdgeEnd) it.next();
    // JTS:       Label label = e.getLabel();
    // JTS:       for (int geomi = 0; geomi < 2; geomi++) {
    // JTS:         if (label.isLine(geomi) && label.getLocation(geomi) == Location.BOUNDARY)
    // JTS:           hasDimensionalCollapseEdge[geomi] = true;
    // JTS:       }
    // JTS:     }
    // JTS: //Debug.print(this);
    // JTS:     for (Iterator it = iterator(); it.hasNext(); ) {
    // JTS:       EdgeEnd e = (EdgeEnd) it.next();
    // JTS:       Label label = e.getLabel();
    // JTS: //Debug.println(e);
    // JTS:       for (int geomi = 0; geomi < 2; geomi++) {
    // JTS:         if (label.isAnyNull(geomi)) {
    // JTS:           int loc = Location.NONE;
    // JTS:           if (hasDimensionalCollapseEdge[geomi]) {
    // JTS:             loc = Location.EXTERIOR;
    // JTS:           }
    // JTS:           else {
    // JTS:             Coordinate p = e.getCoordinate();
    // JTS:             loc = getLocation(geomi, p, geomGraph);
    // JTS:           }
    // JTS:           label.setAllLocationsIfNull(geomi, loc);
    // JTS:         }
    // JTS:       }
    // JTS: //Debug.println(e);
    // JTS:     }
    // JTS: //Debug.print(this);
    // JTS: //Debug.printIfWatch(this);
    // JTS:   }
    pub(crate) fn compute_labeling(&self, graph_a: &GeometryGraph<F>, graph_b: &GeometryGraph<F>) {
        todo!()
    }
}
