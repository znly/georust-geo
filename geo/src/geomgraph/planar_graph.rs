// JTS: import org.locationtech.jts.algorithm.Orientation;
// JTS: import org.locationtech.jts.geom.Coordinate;
// JTS: import org.locationtech.jts.geom.Location;
// JTS: 
// JTS: /**
// JTS:  * The computation of the <code>IntersectionMatrix</code> relies on the use of a structure
// JTS:  * called a "topology graph".  The topology graph contains nodes and edges
// JTS:  * corresponding to the nodes and line segments of a <code>Geometry</code>. Each
// JTS:  * node and edge in the graph is labeled with its topological location relative to
// JTS:  * the source geometry.
// JTS:  * <P>
// JTS:  * Note that there is no requirement that points of self-intersection be a vertex.
// JTS:  * Thus to obtain a correct topology graph, <code>Geometry</code>s must be
// JTS:  * self-noded before constructing their graphs.
// JTS:  * <P>
// JTS:  * Two fundamental operations are supported by topology graphs:
// JTS:  * <UL>
// JTS:  *   <LI>Computing the intersections between all the edges and nodes of a single graph
// JTS:  *   <LI>Computing the intersections between the edges and nodes of two different graphs
// JTS:  * </UL>
// JTS:  *
// JTS:  * @version 1.7
// JTS:  */
// JTS: public class PlanarGraph
// JTS: {
/// The computation of the [IntersectionMatrix] relies on the use of a
/// structure called a "topology graph". The topology graph contains nodes and
/// edges corresponding to the nodes and line segments of a [Geometry]. Each
/// node and edge in the graph is labeled with its topological location
/// relative to the source geometry.
/// 
/// Note that there is no requirement that points of self-intersection be a
/// vertex.  Thus to obtain a correct topology graph, [Geometry] must be
/// self-noded before constructing their graphs.
/// 
/// Two fundamental operations are supported by topology graphs:
///   - Computing the intersections between all the edges and nodes of a single graph
///   - Computing the intersections between the edges and nodes of two different graphs
pub struct PlanarGraph();

impl {
    // JTS:   /**
    // JTS:    * For nodes in the Collection, link the DirectedEdges at the node that are in the result.
    // JTS:    * This allows clients to link only a subset of nodes in the graph, for
    // JTS:    * efficiency (because they know that only a subset is of interest).
    // JTS:    */
    // JTS:   public static void linkResultDirectedEdges(Collection nodes)
    // JTS:   {
    // JTS:     for (Iterator nodeit = nodes.iterator(); nodeit.hasNext(); ) {
    // JTS:       Node node = (Node) nodeit.next();
    // JTS:       ((DirectedEdgeStar) node.getEdges()).linkResultDirectedEdges();
    // JTS:     }
    // JTS:   }
    // JTS: 
    // JTS:   protected List edges        = new ArrayList();
    // JTS:   protected NodeMap nodes;
    // JTS:   protected List edgeEndList  = new ArrayList();
    // JTS: 
    // JTS:   public PlanarGraph(NodeFactory nodeFact) {
    // JTS:     nodes = new NodeMap(nodeFact);
    // JTS:   }
    // JTS: 
    // JTS:   public PlanarGraph() {
    // JTS:     nodes = new NodeMap(new NodeFactory());
    // JTS:   }
    // JTS: 
    // JTS:   public Iterator getEdgeIterator() { return edges.iterator(); }
    // JTS:   public Collection getEdgeEnds() { return edgeEndList; }
    // JTS: 
    // JTS:   public boolean isBoundaryNode(int geomIndex, Coordinate coord)
    // JTS:   {
    // JTS:     Node node = nodes.find(coord);
    // JTS:     if (node == null) return false;
    // JTS:     Label label = node.getLabel();
    // JTS:     if (label != null && label.getLocation(geomIndex) == Location.BOUNDARY) return true;
    // JTS:     return false;
    // JTS:   }
    // JTS:   protected void insertEdge(Edge e)
    // JTS:   {
    // JTS:     edges.add(e);
    // JTS:   }
    // JTS:   public void add(EdgeEnd e)
    // JTS:   {
    // JTS:     nodes.add(e);
    // JTS:     edgeEndList.add(e);
    // JTS:   }
    // JTS: 
    // JTS:   public Iterator getNodeIterator() { return nodes.iterator(); }
    // JTS:   public Collection getNodes() { return nodes.values(); }
    // JTS:   public Node addNode(Node node) { return nodes.addNode(node); }
    // JTS:   public Node addNode(Coordinate coord) { return nodes.addNode(coord); }
    // JTS:   /**
    // JTS:    * @return the node if found; null otherwise
    // JTS:    */
    // JTS:   public Node find(Coordinate coord) { return nodes.find(coord); }
    // JTS: 
    // JTS:   /**
    // JTS:    * Add a set of edges to the graph.  For each edge two DirectedEdges
    // JTS:    * will be created.  DirectedEdges are NOT linked by this method.
    // JTS:    */
    // JTS:   public void addEdges(List edgesToAdd)
    // JTS:   {
    // JTS:     // create all the nodes for the edges
    // JTS:     for (Iterator it = edgesToAdd.iterator(); it.hasNext(); ) {
    // JTS:       Edge e = (Edge) it.next();
    // JTS:       edges.add(e);
    // JTS: 
    // JTS:       DirectedEdge de1 = new DirectedEdge(e, true);
    // JTS:       DirectedEdge de2 = new DirectedEdge(e, false);
    // JTS:       de1.setSym(de2);
    // JTS:       de2.setSym(de1);
    // JTS: 
    // JTS:       add(de1);
    // JTS:       add(de2);
    // JTS:     }
    // JTS:   }
    // JTS: 
    // JTS:   /**
    // JTS:    * Link the DirectedEdges at the nodes of the graph.
    // JTS:    * This allows clients to link only a subset of nodes in the graph, for
    // JTS:    * efficiency (because they know that only a subset is of interest).
    // JTS:    */
    // JTS:   public void linkResultDirectedEdges()
    // JTS:   {
    // JTS:     for (Iterator nodeit = nodes.iterator(); nodeit.hasNext(); ) {
    // JTS:       Node node = (Node) nodeit.next();
    // JTS:       ((DirectedEdgeStar) node.getEdges()).linkResultDirectedEdges();
    // JTS:     }
    // JTS:   }
    // JTS:   /**
    // JTS:    * Link the DirectedEdges at the nodes of the graph.
    // JTS:    * This allows clients to link only a subset of nodes in the graph, for
    // JTS:    * efficiency (because they know that only a subset is of interest).
    // JTS:    */
    // JTS:   public void linkAllDirectedEdges()
    // JTS:   {
    // JTS:     for (Iterator nodeit = nodes.iterator(); nodeit.hasNext(); ) {
    // JTS:       Node node = (Node) nodeit.next();
    // JTS:       ((DirectedEdgeStar) node.getEdges()).linkAllDirectedEdges();
    // JTS:     }
    // JTS:   }
    // JTS:   /**
    // JTS:    * Returns the EdgeEnd which has edge e as its base edge
    // JTS:    * (MD 18 Feb 2002 - this should return a pair of edges)
    // JTS:    *
    // JTS:    * @return the edge, if found
    // JTS:    *    <code>null</code> if the edge was not found
    // JTS:    */
    // JTS:   public EdgeEnd findEdgeEnd(Edge e)
    // JTS:   {
    // JTS:     for (Iterator i = getEdgeEnds().iterator(); i.hasNext(); ) {
    // JTS:       EdgeEnd ee = (EdgeEnd) i.next();
    // JTS:       if (ee.getEdge() == e)
    // JTS:         return ee;
    // JTS:     }
    // JTS:     return null;
    // JTS:   }
    // JTS: 
    // JTS:   /**
    // JTS:    * Returns the edge whose first two coordinates are p0 and p1
    // JTS:    *
    // JTS:    * @return the edge, if found
    // JTS:    *    <code>null</code> if the edge was not found
    // JTS:    */
    // JTS:   public Edge findEdge(Coordinate p0, Coordinate p1)
    // JTS:   {
    // JTS:     for (int i = 0; i < edges.size(); i++) {
    // JTS:       Edge e = (Edge) edges.get(i);
    // JTS:       Coordinate[] eCoord = e.getCoordinates();
    // JTS:       if (p0.equals(eCoord[0]) && p1.equals(eCoord[1]) )
    // JTS:         return e;
    // JTS:     }
    // JTS:     return null;
    // JTS:   }
    // JTS:   /**
    // JTS:    * Returns the edge which starts at p0 and whose first segment is
    // JTS:    * parallel to p1
    // JTS:    *
    // JTS:    * @return the edge, if found
    // JTS:    *    <code>null</code> if the edge was not found
    // JTS:    */
    // JTS:   public Edge findEdgeInSameDirection(Coordinate p0, Coordinate p1)
    // JTS:   {
    // JTS:     for (int i = 0; i < edges.size(); i++) {
    // JTS:       Edge e = (Edge) edges.get(i);
    // JTS: 
    // JTS:       Coordinate[] eCoord = e.getCoordinates();
    // JTS:       if (matchInSameDirection(p0, p1, eCoord[0], eCoord[1]) )
    // JTS:         return e;
    // JTS: 
    // JTS:       if (matchInSameDirection(p0, p1, eCoord[eCoord.length - 1], eCoord[eCoord.length - 2]) )
    // JTS:         return e;
    // JTS:     }
    // JTS:     return null;
    // JTS:   }
    // JTS: 
    // JTS:   /**
    // JTS:    * The coordinate pairs match if they define line segments lying in the same direction.
    // JTS:    * E.g. the segments are parallel and in the same quadrant
    // JTS:    * (as opposed to parallel and opposite!).
    // JTS:    */
    // JTS:   private boolean matchInSameDirection(Coordinate p0, Coordinate p1, Coordinate ep0, Coordinate ep1)
    // JTS:   {
    // JTS:     if (! p0.equals(ep0))
    // JTS:       return false;
    // JTS: 
    // JTS:     if (Orientation.index(p0, p1, ep1) == Orientation.COLLINEAR
    // JTS:          && Quadrant.quadrant(p0, p1) == Quadrant.quadrant(ep0, ep1) )
    // JTS:       return true;
    // JTS:     return false;
    // JTS:   }
    // JTS: 
    // JTS:   public void printEdges(PrintStream out)
    // JTS:   {
    // JTS:     out.println("Edges:");
// JTS:     for (int i = 0; i < edges.size(); i++) {
// JTS:       out.println("edge " + i + ":");
// JTS:       Edge e = (Edge) edges.get(i);
// JTS:       e.print(out);
// JTS:       e.eiList.print(out);
// JTS:     }
// JTS:   }
// JTS:   void debugPrint(Object o)
// JTS:   {
// JTS:     System.out.print(o);
// JTS:   }
// JTS:   void debugPrintln(Object o)
// JTS:   {
// JTS:     System.out.println(o);
// JTS:   }
// JTS: 
// JTS: }
}
