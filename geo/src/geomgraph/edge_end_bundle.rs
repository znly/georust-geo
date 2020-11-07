use super::{EdgeEnd, Float};

// NOTE: in JTS this is in the algorithm::relate package, but because
// we moved EdgeEndBundleStar into geomgraph, we also must move EdgeEndBundle
// (we moved EdgeEndBundleStar rather than wrestle with inheritance of EdgeEndBundle)

// JTS: import org.locationtech.jts.algorithm.BoundaryNodeRule;
// JTS: import org.locationtech.jts.geom.IntersectionMatrix;
// JTS: import org.locationtech.jts.geom.Location;
// JTS: import org.locationtech.jts.geomgraph.Edge;
// JTS: import org.locationtech.jts.geomgraph.EdgeEnd;
// JTS: import org.locationtech.jts.geomgraph.GeometryGraph;
// JTS: import org.locationtech.jts.geomgraph.Label;
// JTS: import org.locationtech.jts.geomgraph.Position;
// JTS:
// JTS: /**
// JTS:  * A collection of {@link EdgeEnd}s which obey the following invariant:
// JTS:  * They originate at the same node and have the same direction.
// JTS:  *
// JTS:  * @version 1.7
// JTS:  */
// JTS: public class EdgeEndBundle
// JTS:   extends EdgeEnd
#[derive(Clone)]
pub(crate) struct EdgeEndBundle<F>
where
    F: Float,
{
    edge_ends: Vec<EdgeEnd<F>>,
}

impl<F> EdgeEndBundle<F>
where
    F: Float,
{
    // JTS: {
    // JTS: //  private BoundaryNodeRule boundaryNodeRule;
    // JTS:   private List edgeEnds = new ArrayList();
    // JTS:
    // JTS:   public EdgeEndBundle(BoundaryNodeRule boundaryNodeRule, EdgeEnd e)
    // JTS:   {
    // JTS:     super(e.getEdge(), e.getCoordinate(), e.getDirectedCoordinate(), new Label(e.getLabel()));
    // JTS:     insert(e);
    // JTS:     /*
    // JTS:     if (boundaryNodeRule != null)
    // JTS:       this.boundaryNodeRule = boundaryNodeRule;
    // JTS:     else
    // JTS:       boundaryNodeRule = BoundaryNodeRule.OGC_SFS_BOUNDARY_RULE;
    // JTS:     */
    // JTS:   }
    // JTS:
    // JTS:   public EdgeEndBundle(EdgeEnd e)
    // JTS:   {
    // JTS:     this(null, e);
    // JTS:   }
    pub(crate) fn new() -> Self {
        // REVIEW: there's a lot more to the JTS initializer since EdgeEndBundle inherits from
        // EdgeEnd, but let's see if it's necessary.
        Self { edge_ends: vec![] }
    }

    // JTS:   public Label getLabel() { return label; }
    // JTS:   public Iterator iterator() { return edgeEnds.iterator(); }
    // JTS:   public List getEdgeEnds() { return edgeEnds; }
    // JTS:
    // JTS:   public void insert(EdgeEnd e)
    // JTS:   {
    // JTS:     // Assert: start point is the same
    // JTS:     // Assert: direction is the same
    // JTS:     edgeEnds.add(e);
    // JTS:   }
    pub(crate) fn insert(&mut self, edge_end: EdgeEnd<F>) {
        self.edge_ends.push(edge_end);
    }

    // JTS:   /**
    // JTS:    * This computes the overall edge label for the set of
    // JTS:    * edges in this EdgeStubBundle.  It essentially merges
    // JTS:    * the ON and side labels for each edge.  These labels must be compatible
    // JTS:    */
    // JTS:   public void computeLabel(BoundaryNodeRule boundaryNodeRule)
    // JTS:   {
    // JTS:     // create the label.  If any of the edges belong to areas,
    // JTS:     // the label must be an area label
    // JTS:     boolean isArea = false;
    // JTS:     for (Iterator it = iterator(); it.hasNext(); ) {
    // JTS:       EdgeEnd e = (EdgeEnd) it.next();
    // JTS:       if (e.getLabel().isArea()) isArea = true;
    // JTS:     }
    // JTS:     if (isArea)
    // JTS:       label = new Label(Location.NONE, Location.NONE, Location.NONE);
    // JTS:     else
    // JTS:       label = new Label(Location.NONE);
    // JTS:
    // JTS:     // compute the On label, and the side labels if present
    // JTS:     for (int i = 0; i < 2; i++) {
    // JTS:       computeLabelOn(i, boundaryNodeRule);
    // JTS:       if (isArea)
    // JTS:         computeLabelSides(i);
    // JTS:     }
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Compute the overall ON location for the list of EdgeStubs.
    // JTS:    * (This is essentially equivalent to computing the self-overlay of a single Geometry)
    // JTS:    * edgeStubs can be either on the boundary (e.g. Polygon edge)
    // JTS:    * OR in the interior (e.g. segment of a LineString)
    // JTS:    * of their parent Geometry.
    // JTS:    * In addition, GeometryCollections use a {@link BoundaryNodeRule} to determine
    // JTS:    * whether a segment is on the boundary or not.
    // JTS:    * Finally, in GeometryCollections it can occur that an edge is both
    // JTS:    * on the boundary and in the interior (e.g. a LineString segment lying on
    // JTS:    * top of a Polygon edge.) In this case the Boundary is given precedence.
    // JTS:    * <br>
    // JTS:    * These observations result in the following rules for computing the ON location:
    // JTS:    * <ul>
    // JTS:    * <li> if there are an odd number of Bdy edges, the attribute is Bdy
    // JTS:    * <li> if there are an even number >= 2 of Bdy edges, the attribute is Int
    // JTS:    * <li> if there are any Int edges, the attribute is Int
    // JTS:    * <li> otherwise, the attribute is NULL.
    // JTS:    * </ul>
    // JTS:    */
    // JTS:   private void computeLabelOn(int geomIndex, BoundaryNodeRule boundaryNodeRule)
    // JTS:   {
    // JTS:     // compute the ON location value
    // JTS:     int boundaryCount = 0;
    // JTS:     boolean foundInterior = false;
    // JTS:
    // JTS:     for (Iterator it = iterator(); it.hasNext(); ) {
    // JTS:       EdgeEnd e = (EdgeEnd) it.next();
    // JTS:       int loc = e.getLabel().getLocation(geomIndex);
    // JTS:       if (loc == Location.BOUNDARY) boundaryCount++;
    // JTS:       if (loc == Location.INTERIOR) foundInterior = true;
    // JTS:     }
    // JTS:     int loc = Location.NONE;
    // JTS:     if (foundInterior)  loc = Location.INTERIOR;
    // JTS:     if (boundaryCount > 0) {
    // JTS:       loc = GeometryGraph.determineBoundary(boundaryNodeRule, boundaryCount);
    // JTS:     }
    // JTS:     label.setLocation(geomIndex, loc);
    // JTS:
    // JTS:   }
    // JTS:   /**
    // JTS:    * Compute the labelling for each side
    // JTS:    */
    // JTS:   private void computeLabelSides(int geomIndex)
    // JTS:   {
    // JTS:     computeLabelSide(geomIndex, Position.LEFT);
    // JTS:     computeLabelSide(geomIndex, Position.RIGHT);
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * To compute the summary label for a side, the algorithm is:
    // JTS:    *   FOR all edges
    // JTS:    *     IF any edge's location is INTERIOR for the side, side location = INTERIOR
    // JTS:    *     ELSE IF there is at least one EXTERIOR attribute, side location = EXTERIOR
    // JTS:    *     ELSE  side location = NULL
    // JTS:    *  <br>
    // JTS:    *  Note that it is possible for two sides to have apparently contradictory information
    // JTS:    *  i.e. one edge side may indicate that it is in the interior of a geometry, while
    // JTS:    *  another edge side may indicate the exterior of the same geometry.  This is
    // JTS:    *  not an incompatibility - GeometryCollections may contain two Polygons that touch
    // JTS:    *  along an edge.  This is the reason for Interior-primacy rule above - it
    // JTS:    *  results in the summary label having the Geometry interior on <b>both</b> sides.
    // JTS:    */
    // JTS:   private void computeLabelSide(int geomIndex, int side)
    // JTS:   {
    // JTS:     for (Iterator it = iterator(); it.hasNext(); ) {
    // JTS:       EdgeEnd e = (EdgeEnd) it.next();
    // JTS:       if (e.getLabel().isArea()) {
    // JTS:         int loc = e.getLabel().getLocation(geomIndex, side);
    // JTS:         if (loc == Location.INTERIOR) {
    // JTS:             label.setLocation(geomIndex, side, Location.INTERIOR);
    // JTS:             return;
    // JTS:         }
    // JTS:         else if (loc == Location.EXTERIOR)
    // JTS:               label.setLocation(geomIndex, side, Location.EXTERIOR);
    // JTS:       }
    // JTS:     }
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Update the IM with the contribution for the computed label for the EdgeStubs.
    // JTS:    */
    // JTS:   void updateIM(IntersectionMatrix im)
    // JTS:   {
    // JTS:     Edge.updateIM(label, im);
    // JTS:   }
    // JTS:   public void print(PrintStream out)
    // JTS:   {
    // JTS:     out.println("EdgeEndBundle--> Label: " + label);
    // JTS:     for (Iterator it = iterator(); it.hasNext(); ) {
    // JTS:       EdgeEnd ee = (EdgeEnd) it.next();
    // JTS:       ee.print(out);
    // JTS:       out.println();
    // JTS:     }
    // JTS:   }
    // JTS: }
}
