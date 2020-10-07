use super::{Coordinate, EdgeIntersection};

use std::collections::BTreeSet;

// JTS: /**
// JTS:  * A list of edge intersections along an {@link Edge}.
// JTS:  * Implements splitting an edge with intersections
// JTS:  * into multiple resultant edges.
// JTS:  *
// JTS:  * @version 1.7
// JTS:  */
// JTS: public class EdgeIntersectionList
// JTS: {
#[derive(PartialEq)]
pub struct EdgeIntersectionList<F: num_traits::Float> {
    node_set: BTreeSet<EdgeIntersection<F>>,
}

impl<F: num_traits::Float> EdgeIntersectionList<F> {
    // JTS:   // a Map <EdgeIntersection, EdgeIntersection>
    // JTS:   private Map nodeMap = new TreeMap();
    // JTS:   Edge edge;  // the parent edge
    // JTS:
    // JTS:   public EdgeIntersectionList(Edge edge)
    // JTS:   {
    // JTS:     this.edge = edge;
    // JTS:   }

    // JTS:   /**
    // JTS:    * Adds an intersection into the list, if it isn't already there.
    // JTS:    * The input segmentIndex and dist are expected to be normalized.
    // JTS:    * @return the EdgeIntersection found or added
    // JTS:    */
    // JTS:   public EdgeIntersection add(Coordinate intPt, int segmentIndex, double dist)
    // JTS:   {
    /// Adds an intersection into the list, if it isn't already there.
    /// The input segmentIndex and dist are expected to be normalized.
    /// @return the EdgeIntersection found or added
    pub fn add(&mut self, intersection_point: Coordinate<F>, segment_index: usize, dist: F) {
        // JTS:     EdgeIntersection eiNew = new EdgeIntersection(intPt, segmentIndex, dist);
        // JTS:     EdgeIntersection ei = (EdgeIntersection) nodeMap.get(eiNew);
        // JTS:     if (ei != null) {
        // JTS:       return ei;
        // JTS:     }
        // JTS:     nodeMap.put(eiNew, eiNew);
        // JTS:     return eiNew;
        let edge_intersection = EdgeIntersection::new(intersection_point, segment_index, dist);
        // BTreeSet only updates the element if it's not alread present
        self.node_set.insert(edge_intersection);

        // Note: the JTS implementation returns the new EdgeIntersection, but it seems unused.
        // Returning it would require some reference gymnastics, so I'm going to omit it until such
        // a time as its needed.

        // JTS:   }
    }
}

// JTS:   /**
// JTS:    * Returns an iterator of {@link EdgeIntersection}s
// JTS:    *
// JTS:    * @return an Iterator of EdgeIntersections
// JTS:    */
// JTS:   public Iterator iterator() { return nodeMap.values().iterator(); }
// JTS:
// JTS:   /**
// JTS:    * Tests if the given point is an edge intersection
// JTS:    *
// JTS:    * @param pt the point to test
// JTS:    * @return true if the point is an intersection
// JTS:    */
// JTS:   public boolean isIntersection(Coordinate pt)
// JTS:   {
// JTS:     for (Iterator it = iterator(); it.hasNext(); ) {
// JTS:       EdgeIntersection ei = (EdgeIntersection) it.next();
// JTS:       if (ei.coord.equals(pt))
// JTS:        return true;
// JTS:     }
// JTS:     return false;
// JTS:   }
// JTS:
// JTS:   /**
// JTS:    * Adds entries for the first and last points of the edge to the list
// JTS:    */
// JTS:   public void addEndpoints()
// JTS:   {
// JTS:     int maxSegIndex = edge.pts.length - 1;
// JTS:     add(edge.pts[0], 0, 0.0);
// JTS:     add(edge.pts[maxSegIndex], maxSegIndex, 0.0);
// JTS:   }
// JTS:
// JTS:   /**
// JTS:    * Creates new edges for all the edges that the intersections in this
// JTS:    * list split the parent edge into.
// JTS:    * Adds the edges to the input list (this is so a single list
// JTS:    * can be used to accumulate all split edges for a Geometry).
// JTS:    *
// JTS:    * @param edgeList a list of EdgeIntersections
// JTS:    */
// JTS:   public void addSplitEdges(List edgeList)
// JTS:   {
// JTS:     // ensure that the list has entries for the first and last point of the edge
// JTS:     addEndpoints();
// JTS:
// JTS:     Iterator it = iterator();
// JTS:     // there should always be at least two entries in the list
// JTS:     EdgeIntersection eiPrev = (EdgeIntersection) it.next();
// JTS:     while (it.hasNext()) {
// JTS:       EdgeIntersection ei = (EdgeIntersection) it.next();
// JTS:       Edge newEdge = createSplitEdge(eiPrev, ei);
// JTS:       edgeList.add(newEdge);
// JTS:
// JTS:       eiPrev = ei;
// JTS:     }
// JTS:   }
// JTS:   /**
// JTS:    * Create a new "split edge" with the section of points between
// JTS:    * (and including) the two intersections.
// JTS:    * The label for the new edge is the same as the label for the parent edge.
// JTS:    */
// JTS:   Edge createSplitEdge(EdgeIntersection ei0, EdgeIntersection ei1)
// JTS:   {
// JTS: //Debug.print("\ncreateSplitEdge"); Debug.print(ei0); Debug.print(ei1);
// JTS:     int npts = ei1.segmentIndex - ei0.segmentIndex + 2;
// JTS:
// JTS:     Coordinate lastSegStartPt = edge.pts[ei1.segmentIndex];
// JTS:     // if the last intersection point is not equal to the its segment start pt,
// JTS:     // add it to the points list as well.
// JTS:     // (This check is needed because the distance metric is not totally reliable!)
// JTS:     // The check for point equality is 2D only - Z values are ignored
// JTS:     boolean useIntPt1 = ei1.dist > 0.0 || ! ei1.coord.equals2D(lastSegStartPt);
// JTS:     if (! useIntPt1) {
// JTS:       npts--;
// JTS:     }
// JTS:
// JTS:     Coordinate[] pts = new Coordinate[npts];
// JTS:     int ipt = 0;
// JTS:     pts[ipt++] = new Coordinate(ei0.coord);
// JTS:     for (int i = ei0.segmentIndex + 1; i <= ei1.segmentIndex; i++) {
// JTS:       pts[ipt++] = edge.pts[i];
// JTS:     }
// JTS:     if (useIntPt1) pts[ipt] = ei1.coord;
// JTS:     return new Edge(pts, new Label(edge.label));
// JTS:   }
// JTS:
// JTS:   public void print(PrintStream out)
// JTS:   {
// JTS:     out.println("Intersections:");
// JTS:     for (Iterator it = iterator(); it.hasNext(); ) {
// JTS:       EdgeIntersection ei = (EdgeIntersection) it.next();
// JTS:       ei.print(out);
// JTS:     }
// JTS:   }
// JTS: }
