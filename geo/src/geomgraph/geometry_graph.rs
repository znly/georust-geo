// JTS: import java.util.HashMap;
// JTS: import java.util.Iterator;
// JTS: import java.util.List;
// JTS: import java.util.Map;

use super::{
    algorithm::{
        boundary_node_rule::{BoundaryNodeRule, Mod2BoundaryNodeRule},
        LineIntersector,
    },
    index::{EdgeSetIntersector, SegmentIntersector, SimpleEdgeSetIntersector},
    BasicNodeFactory, Coordinate, Edge, Float, GraphComponent, Label, Location, Node, NodeFactory,
    PlanarGraph, Position,
};

use crate::algorithm::dimensions::HasDimensions;
use crate::{Geometry, Line, LineString, Point, Polygon};

use std::cell::RefCell;

// JTS: import org.locationtech.jts.algorithm.BoundaryNodeRule;
// JTS: import org.locationtech.jts.algorithm.LineIntersector;
// JTS: import org.locationtech.jts.algorithm.Orientation;
// JTS: import org.locationtech.jts.algorithm.PointLocator;
// JTS: import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
// JTS: import org.locationtech.jts.algorithm.locate.PointOnGeometryLocator;
// JTS: import org.locationtech.jts.geom.Coordinate;
// JTS: import org.locationtech.jts.geom.CoordinateArrays;
// JTS: import org.locationtech.jts.geom.Geometry;
// JTS: import org.locationtech.jts.geom.GeometryCollection;
// JTS: import org.locationtech.jts.geom.LineString;
// JTS: import org.locationtech.jts.geom.LinearRing;
// JTS: import org.locationtech.jts.geom.Location;
// JTS: import org.locationtech.jts.geom.MultiLineString;
// JTS: import org.locationtech.jts.geom.MultiPoint;
// JTS: import org.locationtech.jts.geom.MultiPolygon;
// JTS: import org.locationtech.jts.geom.Point;
// JTS: import org.locationtech.jts.geom.Polygon;
// JTS: import org.locationtech.jts.geom.Polygonal;
// JTS: import org.locationtech.jts.geomgraph.index.EdgeSetIntersector;
// JTS: import org.locationtech.jts.geomgraph.index.SegmentIntersector;
// JTS: import org.locationtech.jts.geomgraph.index.SimpleMCSweepLineIntersector;
// JTS: import org.locationtech.jts.util.Assert;

// JTS: /**
// JTS:  * A GeometryGraph is a graph that models a given Geometry
// JTS:  * @version 1.7
// JTS:  */
// JTS: public class GeometryGraph
// JTS:   extends PlanarGraph
// JTS: {
/// A GeometryGraph is a graph that models a given Geometry
pub(crate) struct GeometryGraph<'a, F>
where
    F: Float,
{
    arg_index: usize,
    parent_geometry: &'a Geometry<F>,
    use_boundary_determination_rule: bool,
    planar_graph: PlanarGraph<F>,
}

// PlanarGraph delegations
// In JTS this is achieved through inheritance - GeometryGraph inherits from PlanarGraph
impl<F> GeometryGraph<'_, F>
where
    F: Float,
{
    pub fn edges(&self) -> &[RefCell<Edge<F>>] {
        self.planar_graph.edges()
    }

    pub fn insert_edge(&mut self, edge: Edge<F>) {
        self.planar_graph.insert_edge(edge)
    }

    pub fn is_boundary_node(&self, geom_index: usize, coord: Coordinate<F>) -> bool {
        self.planar_graph.is_boundary_node(geom_index, coord)
    }

    pub fn add_node_with_coordinate(&mut self, coord: Coordinate<F>) -> &mut Node<F> {
        self.planar_graph.add_node_with_coordinate(coord)
    }

    pub fn nodes_iter(&self) -> impl Iterator<Item = &Node<F>> {
        self.planar_graph.nodes.iter()
    }
}

impl<'a, F> GeometryGraph<'a, F>
where
    F: Float,
{
    // JTS: /**
    // JTS:  * This method implements the Boundary Determination Rule
    // JTS:  * for determining whether
    // JTS:  * a component (node or edge) that appears multiple times in elements
    // JTS:  * of a MultiGeometry is in the boundary or the interior of the Geometry
    // JTS:  * <br>
    // JTS:  * The SFS uses the "Mod-2 Rule", which this function implements
    // JTS:  * <br>
    // JTS:  * An alternative (and possibly more intuitive) rule would be
    // JTS:  * the "At Most One Rule":
    // JTS:  *    isInBoundary = (componentCount == 1)
    // JTS:  */
    // JTS: /*
    // JTS:   public static boolean isInBoundary(int boundaryCount)
    // JTS:   {
    // JTS:     // the "Mod-2 Rule"
    // JTS:     return boundaryCount % 2 == 1;
    // JTS:   }
    // JTS:   public static int determineBoundary(int boundaryCount)
    // JTS:   {
    // JTS:     return isInBoundary(boundaryCount) ? Location.BOUNDARY : Location.INTERIOR;
    // JTS:   }
    // JTS: */
    // JTS:
    // JTS:   public static int determineBoundary(BoundaryNodeRule boundaryNodeRule, int boundaryCount)
    // JTS:   {
    // JTS:     return boundaryNodeRule.isInBoundary(boundaryCount)
    // JTS:         ? Location.BOUNDARY : Location.INTERIOR;
    // JTS:   }
    pub fn determine_boundary(
        boundary_node_rule: &impl BoundaryNodeRule,
        boundary_count: usize,
    ) -> Location {
        if boundary_node_rule.is_in_boundary(boundary_count) {
            Location::Boundary
        } else {
            Location::Interior
        }
    }

    // JTS:
    // JTS:   private Geometry parentGeom;
    // JTS:
    // JTS:   /**
    // JTS:    * The lineEdgeMap is a map of the linestring components of the
    // JTS:    * parentGeometry to the edges which are derived from them.
    // JTS:    * This is used to efficiently perform findEdge queries
    // JTS:    */
    // JTS:   private Map lineEdgeMap = new HashMap();
    // JTS:
    // JTS:   private BoundaryNodeRule boundaryNodeRule = null;
    // JTS:
    // JTS:   /**
    // JTS:    * If this flag is true, the Boundary Determination Rule will used when deciding
    // JTS:    * whether nodes are in the boundary or not
    // JTS:    */
    // JTS:   private boolean useBoundaryDeterminationRule = true;
    // JTS:   private int argIndex;  // the index of this geometry as an argument to a spatial function (used for labelling)
    // JTS:   private Collection boundaryNodes;
    // JTS:   private boolean hasTooFewPoints = false;
    // JTS:   private Coordinate invalidPoint = null;
    // JTS:
    // JTS:   private PointOnGeometryLocator areaPtLocator = null;
    // JTS:   // for use if geometry is not Polygonal
    // JTS:   private final PointLocator ptLocator = new PointLocator();
    // JTS:
    // JTS:   private EdgeSetIntersector createEdgeSetIntersector()
    // JTS:   {
    // CLEANUP: hardcode the type? is perf significant?
    fn create_edge_set_intersector() -> Box<dyn EdgeSetIntersector<F>> {
        // JTS:   // various options for computing intersections, from slowest to fastest
        // JTS:
        // JTS:   //private EdgeSetIntersector esi = new SimpleEdgeSetIntersector();
        // JTS:   //private EdgeSetIntersector esi = new MonotoneChainIntersector();
        // JTS:   //private EdgeSetIntersector esi = new NonReversingChainIntersector();
        // JTS:   //private EdgeSetIntersector esi = new SimpleSweepLineIntersector();
        // JTS:   //private EdgeSetIntersector esi = new MCSweepLineIntersector();
        // JTS:
        // JTS:     //return new SimpleEdgeSetIntersector();
        // TODO: implement more optimized edge set intersectors
        return Box::new(SimpleEdgeSetIntersector::new());
        // JTS:     return new SimpleMCSweepLineIntersector();
        // JTS:   }
    }

    // JTS:   public GeometryGraph(int argIndex, Geometry parentGeom)
    // JTS:   {
    // JTS:     this(argIndex, parentGeom,
    // JTS:          BoundaryNodeRule.OGC_SFS_BOUNDARY_RULE
    // JTS:          );
    // JTS:   }
    // JTS:
    // JTS:   public GeometryGraph(int argIndex, Geometry parentGeom, BoundaryNodeRule boundaryNodeRule) {
    // JTS:     this.argIndex = argIndex;
    // JTS:     this.parentGeom = parentGeom;
    // JTS:     this.boundaryNodeRule = boundaryNodeRule;
    // JTS:     if (parentGeom != null) {
    // JTS: //      precisionModel = parentGeom.getPrecisionModel();
    // JTS: //      SRID = parentGeom.getSRID();
    // JTS:       add(parentGeom);
    // JTS:     }
    // JTS:   }
    pub fn new(arg_index: usize, parent_geometry: &'a Geometry<F>) -> Self {
        let mut graph = GeometryGraph {
            arg_index,
            parent_geometry,
            use_boundary_determination_rule: true,
            planar_graph: PlanarGraph::new(),
        };
        graph.add_geometry(parent_geometry);
        graph
    }

    // JTS:
    // JTS:   /**
    // JTS:    * This constructor is used by clients that wish to add Edges explicitly,
    // JTS:    * rather than adding a Geometry.  (An example is BufferOp).
    // JTS:    */
    // JTS:   // no longer used
    // JTS: //  public GeometryGraph(int argIndex, PrecisionModel precisionModel, int SRID) {
    // JTS: //    this(argIndex, null);
    // JTS: //    this.precisionModel = precisionModel;
    // JTS: //    this.SRID = SRID;
    // JTS: //  }
    // JTS: //  public PrecisionModel getPrecisionModel()
    // JTS: //  {
    // JTS: //    return precisionModel;
    // JTS: //  }
    // JTS: //  public int getSRID() { return SRID; }
    // JTS:
    // JTS:   public boolean hasTooFewPoints() { return hasTooFewPoints; }
    // JTS:
    // JTS:   public Coordinate getInvalidPoint() { return invalidPoint; }
    // JTS:
    // JTS:   public Geometry getGeometry() { return parentGeom; }
    pub fn geometry(&self) -> &Geometry<F> {
        self.parent_geometry
    }

    // JTS:   public BoundaryNodeRule getBoundaryNodeRule() { return boundaryNodeRule; }
    // JTS:
    // JTS:   public Collection getBoundaryNodes()
    // JTS:   {
    // JTS:     if (boundaryNodes == null)
    // JTS:       boundaryNodes = nodes.getBoundaryNodes(argIndex);
    // JTS:     return boundaryNodes;
    // JTS:   }
    fn boundary_nodes(&self) -> Vec<&Node<F>> {
        // TODO: should we need to memoize this like JTS does?
        self.planar_graph.nodes.boundary_nodes(self.arg_index)
    }
    // JTS:
    // JTS:   public Coordinate[] getBoundaryPoints()
    // JTS:   {
    // JTS:     Collection coll = getBoundaryNodes();
    // JTS:     Coordinate[] pts = new Coordinate[coll.size()];
    // JTS:     int i = 0;
    // JTS:     for (Iterator it = coll.iterator(); it.hasNext(); ) {
    // JTS:       Node node = (Node) it.next();
    // JTS:       pts[i++] = node.getCoordinate().copy();
    // JTS:     }
    // JTS:     return pts;
    // JTS:   }
    // JTS:
    // JTS:   public Edge findEdge(LineString line)
    // JTS:   {
    // JTS:     return (Edge) lineEdgeMap.get(line);
    // JTS:   }
    // JTS:
    // JTS:   public void computeSplitEdges(List edgelist)
    // JTS:   {
    // JTS:     for (Iterator i = edges.iterator(); i.hasNext(); ) {
    // JTS:       Edge e = (Edge) i.next();
    // JTS:       e.eiList.addSplitEdges(edgelist);
    // JTS:     }
    // JTS:   }

    // JTS:   private void add(Geometry g)
    // JTS:   {
    pub fn add_geometry(&mut self, geometry: &Geometry<F>) {
        // JTS:     if (g.isEmpty()) return;
        if geometry.is_empty() {
            return;
        }
        // JTS:
        // JTS:     if (g instanceof Polygon)                 addPolygon((Polygon) g);
        // JTS:                         // LineString also handles LinearRings
        // JTS:     else if (g instanceof LineString)         addLineString((LineString) g);
        // JTS:     else if (g instanceof Point)              addPoint((Point) g);
        // JTS:     else if (g instanceof MultiPoint)         addCollection((MultiPoint) g);
        // JTS:     else if (g instanceof MultiLineString)    addCollection((MultiLineString) g);
        // JTS:     else if (g instanceof MultiPolygon)       addCollection((MultiPolygon) g);
        // JTS:     else if (g instanceof GeometryCollection) addCollection((GeometryCollection) g);
        // JTS:     else  throw new UnsupportedOperationException(g.getClass().getName());
        // JTS:   }
        match geometry {
            Geometry::Line(line) => self.add_line(line),
            Geometry::Rect(rect) => {
                // PERF: avoid this conversion/clone?
                self.add_polygon(&Polygon::from(rect.clone()));
            }
            Geometry::Point(point) => {
                self.add_point(point);
            }
            Geometry::Polygon(polygon) => self.add_polygon(polygon),
            Geometry::Triangle(triangle) => {
                // PERF: avoid this conversion/clone?
                self.add_polygon(&Polygon::from(triangle.clone()));
            }
            Geometry::LineString(line_string) => self.add_line_string(line_string),
            Geometry::MultiPoint(multi_point) => {
                for point in &multi_point.0 {
                    self.add_point(point);
                }
            }
            Geometry::MultiPolygon(multi_polygon) => {
                // JTS:     // check if this Geometry should obey the Boundary Determination Rule
                // JTS:     // all collections except MultiPolygons obey the rule
                // JTS:     if (g instanceof MultiPolygon)
                // JTS:       useBoundaryDeterminationRule = false;
                self.use_boundary_determination_rule = false;
                for polygon in &multi_polygon.0 {
                    self.add_polygon(polygon);
                }
            }
            Geometry::MultiLineString(multi_line_string) => {
                for line_string in &multi_line_string.0 {
                    // PERF: can we get rid of these clones?
                    self.add_line_string(line_string);
                }
            }
            Geometry::GeometryCollection(geometry_collection) => {
                for geometry in geometry_collection {
                    self.add_geometry(geometry);
                }
            }
        }
    }

    // JTS:   private void addCollection(GeometryCollection gc)
    // JTS:   {
    // JTS:     for (int i = 0; i < gc.getNumGeometries(); i++) {
    // JTS:       Geometry g = gc.getGeometryN(i);
    // JTS:       add(g);
    // JTS:     }
    // JTS:   }
    // JTS:   /**
    // JTS:    * Add a Point to the graph.
    // JTS:    */
    // JTS:   private void addPoint(Point p)
    // JTS:   {
    // JTS:     Coordinate coord = p.getCoordinate();
    // JTS:     insertPoint(argIndex, coord, Location.INTERIOR);
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Adds a polygon ring to the graph.
    // JTS:    * Empty rings are ignored.
    // JTS:    *
    // JTS:    * The left and right topological location arguments assume that the ring is oriented CW.
    // JTS:    * If the ring is in the opposite orientation,
    // JTS:    * the left and right locations must be interchanged.
    // JTS:    */
    // JTS:   private void addPolygonRing(LinearRing lr, int cwLeft, int cwRight)
    // JTS:   {
    fn add_polygon_ring(
        &mut self,
        linear_ring: &LineString<F>,
        cw_left: Location,
        cw_right: Location,
    ) {
        debug_assert!(linear_ring.is_closed());
        // JTS:   	// don't bother adding empty holes
        // JTS:   	if (lr.isEmpty()) return;
        if linear_ring.is_empty() {
            return;
        }

        // JTS:     Coordinate[] coord = CoordinateArrays.removeRepeatedPoints(lr.getCoordinates());
        let mut coords: Vec<Coordinate<F>> = vec![];
        for coord in &linear_ring.0 {
            if coords.last() != Some(coord) {
                coords.push(*coord)
            }
        }

        // JTS:     if (coord.length < 4) {
        // JTS:       hasTooFewPoints = true;
        // JTS:       invalidPoint = coord[0];
        // JTS:       return;
        // JTS:     }
        if coords.len() < 4 {
            todo!("handle invalid ring")
        }
        let first_point = coords[0].clone();

        // JTS:
        // JTS:     int left  = cwLeft;
        // JTS:     int right = cwRight;
        // JTS:     if (Orientation.isCCW(coord)) {
        // JTS:       left = cwRight;
        // JTS:       right = cwLeft;
        // JTS:     }
        use crate::algorithm::winding_order::{Winding, WindingOrder};
        let (left, right) = if linear_ring.winding_order() == Some(WindingOrder::CounterClockwise) {
            (cw_right, cw_left)
        } else {
            (cw_left, cw_right)
        };

        // JTS:     Edge e = new Edge(coord,
        // JTS:                         new Label(argIndex, Location.BOUNDARY, left, right));
        let edge = Edge::new(
            coords,
            Label::new_with_geom_locations(self.arg_index, Location::Boundary, left, right),
        );
        // JTS:     lineEdgeMap.put(lr, e);
        // REVIEW: note we don't implement lineEdgeMap. I don't think we *need* it for the Relate
        // operations and I think it'll require moving edges into a RefCell.

        // JTS:     insertEdge(e);
        self.insert_edge(edge);

        // JTS:     // insert the endpoint as a node, to mark that it is on the boundary
        // JTS:     insertPoint(argIndex, coord[0], Location.BOUNDARY);
        // insert the endpoint as a node, to mark that it is on the boundary
        self.insert_point(self.arg_index, first_point, Location::Boundary);
        // JTS:   }
    }

    // JTS:   private void addPolygon(Polygon p)
    // JTS:   {
    fn add_polygon(&mut self, polygon: &Polygon<F>) {
        // JTS:     addPolygonRing(
        // JTS:             p.getExteriorRing(),
        // JTS:             Location.EXTERIOR,
        // JTS:             Location.INTERIOR);
        self.add_polygon_ring(polygon.exterior(), Location::Exterior, Location::Interior);
        // JTS:     for (int i = 0; i < p.getNumInteriorRing(); i++) {
        // JTS:     	LinearRing hole = p.getInteriorRingN(i);
        // JTS:
        // JTS:       // Holes are topologically labelled opposite to the shell, since
        // JTS:       // the interior of the polygon lies on their opposite side
        // JTS:       // (on the left, if the hole is oriented CW)
        // JTS:       addPolygonRing(
        // JTS:       		hole,
        // JTS:           Location.INTERIOR,
        // JTS:           Location.EXTERIOR);
        // JTS:     }
        // JTS:   }
        for hole in polygon.interiors() {
            self.add_polygon_ring(hole, Location::Interior, Location::Exterior)
        }
    }

    // JTS:   private void addLineString(LineString line)
    // JTS:   {
    fn add_line_string(&mut self, line_string: &LineString<F>) {
        // JTS:     Coordinate[] coord = CoordinateArrays.removeRepeatedPoints(line.getCoordinates());
        let mut coords: Vec<Coordinate<F>> = vec![];
        for coord in &line_string.0 {
            if coords.last() != Some(coord) {
                coords.push(*coord)
            }
        }

        // JTS:     if (coord.length < 2) {
        // JTS:       hasTooFewPoints = true;
        // JTS:       invalidPoint = coord[0];
        // JTS:       return;
        // JTS:     }
        if coords.len() < 2 {
            todo!("handle invalid line string");
        }
        self.insert_boundary_point(self.arg_index, *coords.first().unwrap());
        self.insert_boundary_point(self.arg_index, *coords.last().unwrap());

        // JTS:
        // JTS:     // add the edge for the LineString
        // JTS:     // line edges do not have locations for their left and right sides
        // JTS:     Edge e = new Edge(coord, new Label(argIndex, Location.INTERIOR));
        let edge = Edge::new(
            coords,
            Label::new_with_geom_on_location(self.arg_index, Some(Location::Interior)),
        );

        // REVIEW: note we don't implement lineEdgeMap. I don't think we *need* it for the Relate
        // operations and I think it'll require moving edges into a RefCell.
        // JTS:     lineEdgeMap.put(line, e);

        // JTS:     insertEdge(e);
        self.insert_edge(edge);

        // JTS:     /**
        // JTS:      * Add the boundary points of the LineString, if any.
        // JTS:      * Even if the LineString is closed, add both points as if they were endpoints.
        // JTS:      * This allows for the case that the node already exists and is a boundary point.
        // JTS:      */
        // JTS:     Assert.isTrue(coord.length >= 2, "found LineString with single point");
        // JTS:     insertBoundaryPoint(argIndex, coord[0]);
        // JTS:     insertBoundaryPoint(argIndex, coord[coord.length - 1]);
        // REVIEW: re-ordered code to insert boundary points *before* `coords` moves into Edge::new
        // JTS:   }
    }

    fn add_line(&mut self, line: &Line<F>) {
        self.insert_boundary_point(self.arg_index, line.start);
        self.insert_boundary_point(self.arg_index, line.end);

        let edge = Edge::new(
            vec![line.start, line.end],
            Label::new_with_geom_on_location(self.arg_index, Some(Location::Interior)),
        );

        self.insert_edge(edge);
    }

    // JTS:
    // JTS:   /**
    // JTS:    * Add an Edge computed externally.  The label on the Edge is assumed
    // JTS:    * to be correct.
    // JTS:    */
    // JTS:   public void addEdge(Edge e)
    // JTS:   {
    // JTS:     insertEdge(e);
    // JTS:     Coordinate[] coord = e.getCoordinates();
    // JTS:     // insert the endpoint as a node, to mark that it is on the boundary
    // JTS:     insertPoint(argIndex, coord[0], Location.BOUNDARY);
    // JTS:     insertPoint(argIndex, coord[coord.length - 1], Location.BOUNDARY);
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Add a point computed externally.  The point is assumed to be a
    // JTS:    * Point Geometry part, which has a location of INTERIOR.
    // JTS:    */
    // JTS:   public void addPoint(Coordinate pt)
    // JTS:   {
    // JTS:     insertPoint(argIndex, pt, Location.INTERIOR);
    // JTS:   }
    /// Add a point computed externally.  The point is assumed to be a
    /// Point Geometry part, which has a location of INTERIOR.
    fn add_point(&mut self, point: &Point<F>) {
        self.insert_point(self.arg_index, point.clone().into(), Location::Interior);
    }

    // JTS:   /**
    // JTS:    * Compute self-nodes, taking advantage of the Geometry type to
    // JTS:    * minimize the number of intersection tests.  (E.g. rings are
    // JTS:    * not tested for self-intersection, since they are assumed to be valid).
    // JTS:    *
    // JTS:    * @param li the LineIntersector to use
    // JTS:    * @param computeRingSelfNodes if <code>false</code>, intersection checks are optimized to not test rings for self-intersection
    // JTS:    * @return the computed SegmentIntersector containing information about the intersections found
    // JTS:    */
    // JTS:   public SegmentIntersector computeSelfNodes(LineIntersector li, boolean computeRingSelfNodes)
    // JTS:   {
    // JTS: 	  return computeSelfNodes(li, computeRingSelfNodes, false);
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Compute self-nodes, taking advantage of the Geometry type to
    // JTS:    * minimize the number of intersection tests.  (E.g. rings are
    // JTS:    * not tested for self-intersection, since they are assumed to be valid).
    // JTS:    *
    // JTS:    * @param li the LineIntersector to use
    // JTS:    * @param computeRingSelfNodes if <code>false</code>, intersection checks are optimized to not test rings for self-intersection
    // JTS:    * @param isDoneIfProperInt short-circuit the intersection computation if a proper intersection is found
    // JTS:    * @return the computed SegmentIntersector containing information about the intersections found
    // JTS:    */
    // JTS:   public SegmentIntersector computeSelfNodes(LineIntersector li, boolean computeRingSelfNodes, boolean isDoneIfProperInt)
    // JTS:   {
    /// Compute self-nodes, taking advantage of the Geometry type to minimize the number of
    /// intersection tests.  (E.g. rings are not tested for self-intersection, since they are
    /// assumed to be valid).
    ///
    /// @param li the LineIntersector to use
    /// @param computeRingSelfNodes if <code>false</code>, intersection checks are optimized to not test rings for self-intersection
    /// @param isDoneIfProperInt short-circuit the intersection computation if a proper intersection is found
    /// @return the computed SegmentIntersector containing information about the intersections found
    pub fn compute_self_nodes(
        &mut self,
        line_intersector: Box<dyn LineIntersector<F>>,
        _compute_ring_self_nodes: bool,
        is_done_when_proper_intersection: bool,
    ) -> SegmentIntersector<F> {
        // JTS:     SegmentIntersector si = new SegmentIntersector(li, true, false);
        let mut segment_intersector = SegmentIntersector::new(line_intersector, true, false);

        // JTS:     si.setIsDoneIfProperInt(isDoneIfProperInt);
        segment_intersector.set_is_done_when_proper_intersection(is_done_when_proper_intersection);

        // JTS:     EdgeSetIntersector esi = createEdgeSetIntersector();
        let mut edge_set_intersector = Self::create_edge_set_intersector();

        // TODO: optimize intersection search for valid Polygons and LinearRings
        // JTS:     // optimize intersection search for valid Polygons and LinearRings
        // JTS:     boolean isRings = parentGeom instanceof LinearRing
        // JTS: 			|| parentGeom instanceof Polygon
        // JTS: 			|| parentGeom instanceof MultiPolygon;
        // JTS:     boolean computeAllSegments = computeRingSelfNodes || ! isRings;
        let compute_all_segments = true;

        // JTS:     esi.computeIntersections(edges, si, computeAllSegments);
        edge_set_intersector.compute_intersections(
            self.edges(),
            &mut segment_intersector,
            compute_all_segments,
        );

        // JTS:     //System.out.println("SegmentIntersector # tests = " + si.numTests);
        // JTS:     addSelfIntersectionNodes(argIndex);
        // CLEANUP: read self.arg_index as property within addSelfIntersectionNodes rather than
        //          pass as param?
        self.add_self_intersection_nodes(self.arg_index);

        // JTS:     return si;
        // JTS:   }
        segment_intersector
    }

    // JTS:   public SegmentIntersector computeEdgeIntersections(
    // JTS:     GeometryGraph g,
    // JTS:     LineIntersector li,
    // JTS:     boolean includeProper)
    // JTS:   {
    pub fn compute_edge_intersections(
        &self,
        other: &GeometryGraph<F>,
        line_intersector: Box<dyn LineIntersector<F>>,
        include_proper: bool,
    ) -> SegmentIntersector<F> {
        // JTS:     SegmentIntersector si = new SegmentIntersector(li, includeProper, true);
        // JTS:     si.setBoundaryNodes(this.getBoundaryNodes(), g.getBoundaryNodes());
        let mut segment_intersector =
            SegmentIntersector::new(line_intersector, include_proper, true);
        segment_intersector.set_boundary_nodes(
            // CLEANUP: surely there's a nicer way to: `Vec<&Node> -> Vec<Node>`
            self.boundary_nodes()
                .into_iter()
                .map(|n| n.clone())
                .collect(),
            other
                .boundary_nodes()
                .into_iter()
                .map(|n| n.clone())
                .collect(),
        );

        // JTS:
        // JTS:     EdgeSetIntersector esi = createEdgeSetIntersector();
        // JTS:     esi.computeIntersections(edges, g.edges, si);
        let mut edge_set_intersector = Self::create_edge_set_intersector();
        edge_set_intersector.compute_intersections_testing_all_segments(
            self.edges(),
            other.edges(),
            &mut segment_intersector,
        );

        // JTS: /*
        // JTS: for (Iterator i = g.edges.iterator(); i.hasNext();) {
        // JTS: Edge e = (Edge) i.next();
        // JTS: Debug.print(e.getEdgeIntersectionList());
        // JTS: }
        // JTS: */
        // JTS:     return si;
        segment_intersector
        // JTS:   }
    }

    // JTS:   private void insertPoint(int argIndex, Coordinate coord, int onLocation)
    // JTS:   {
    // JTS:     Node n = nodes.addNode(coord);
    // JTS:     Label lbl = n.getLabel();
    // JTS:     if (lbl == null) {
    // JTS:       n.label = new Label(argIndex, onLocation);
    // JTS:     }
    // JTS:     else
    // JTS:       lbl.setLocation(argIndex, onLocation);
    // JTS:   }
    fn insert_point(&mut self, arg_index: usize, coord: Coordinate<F>, location: Location) {
        let node: &mut Node<F> = self.add_node_with_coordinate(coord);
        // CLEANUP: can we get rid of the Option? Or do we need Edges to maintain Option and share GraphComponent trait
        // VERIFY: JTS does a null check here, but not for boundary points. Can we safely skip it?
        let label: &mut Label = node.label_mut().unwrap();
        label.set_on_location(arg_index, Some(location))
    }

    // JTS:   /**
    // JTS:    * Adds candidate boundary points using the current {@link BoundaryNodeRule}.
    // JTS:    * This is used to add the boundary
    // JTS:    * points of dim-1 geometries (Curves/MultiCurves).
    // JTS:    */
    // JTS:   private void insertBoundaryPoint(int argIndex, Coordinate coord)
    // JTS:   {
    /// Adds candidate boundary points using the current {@link BoundaryNodeRule}.
    /// This is used to add the boundary points of dim-1 geometries (Curves/MultiCurves).
    fn insert_boundary_point(&mut self, arg_index: usize, coord: Coordinate<F>) {
        // JTS:     Node n = nodes.addNode(coord);
        let node: &mut Node<F> = self.add_node_with_coordinate(coord);

        // JTS:     // nodes always have labels
        // JTS:     Label lbl = n.getLabel();
        // nodes always have labels
        // CLEANUP: can we get rid of the Option? Or do we need Edges to maintain Option and share GraphComponent trait
        let label: &mut Label = node.label_mut().unwrap();

        // JTS:     // the new point to insert is on a boundary
        // JTS:     int boundaryCount = 1;
        // JTS:     // determine the current location for the point (if any)
        // JTS:     int loc = Location.NONE;
        // JTS:     loc = lbl.getLocation(argIndex, Position.ON);
        // JTS:     if (loc == Location.BOUNDARY) boundaryCount++;
        // the new point to insert is on a boundary
        let mut boundary_count = 1;
        // determine the current location for the point (if any)
        let location = label.location(arg_index, Position::On);
        if let Some(Location::Boundary) = location {
            boundary_count += 1;
        }

        // JTS:     // determine the boundary status of the point according to the Boundary Determination Rule
        // JTS:     int newLoc = determineBoundary(boundaryNodeRule, boundaryCount);
        // JTS:     lbl.setLocation(argIndex, newLoc);
        // determine the boundary status of the point according to the Boundary Determination Rule
        // TODO: accommodate pluggable boundary node rules
        let new_location = Self::determine_boundary(&Mod2BoundaryNodeRule, boundary_count);
        label.set_on_location(arg_index, Some(new_location));
        // JTS:   }
    }
    // JTS:
    // JTS:   private void addSelfIntersectionNodes(int argIndex)
    // JTS:   {
    fn add_self_intersection_nodes(&mut self, arg_index: usize) {
        // JTS:     for (Iterator i = edges.iterator(); i.hasNext(); ) {
        // JTS:       Edge e = (Edge) i.next();
        // JTS:       int eLoc = e.getLabel().getLocation(argIndex);
        // JTS:       for (Iterator eiIt = e.eiList.iterator(); eiIt.hasNext(); ) {
        // JTS:         EdgeIntersection ei = (EdgeIntersection) eiIt.next();
        // JTS:         addSelfIntersectionNode(argIndex, ei.coord, eLoc);
        // JTS:       }
        // JTS:     }
        // JTS:   }

        let locations_and_intersections: Vec<(Option<Location>, Vec<Coordinate<F>>)> = self
            .edges()
            .into_iter()
            .map(RefCell::borrow)
            .map(|edge| {
                let location = edge.label().and_then(|label| label.on_location(arg_index));

                let coordinates = edge
                    .edge_intersections()
                    .into_iter()
                    .map(|edge_intersection| edge_intersection.coordinate());

                (location, coordinates.collect())
            })
            .collect();

        for (location, edge_intersection_coordinates) in locations_and_intersections {
            for coordinate in edge_intersection_coordinates {
                self.add_self_intersection_node(arg_index, coordinate, location)
            }
        }
    }

    // JTS:   /**
    // JTS:    * Add a node for a self-intersection.
    // JTS:    * If the node is a potential boundary node (e.g. came from an edge which
    // JTS:    * is a boundary) then insert it as a potential boundary node.
    // JTS:    * Otherwise, just add it as a regular node.
    // JTS:    */
    // JTS:   private void addSelfIntersectionNode(int argIndex, Coordinate coord, int loc)
    // JTS:   {
    /// Add a node for a self-intersection.
    ///
    /// If the node is a potential boundary node (e.g. came from an edge which is a boundary), then
    /// insert it as a potential boundary node.  Otherwise, just add it as a regular node.
    fn add_self_intersection_node(
        &mut self,
        arg_index: usize,
        coord: Coordinate<F>,
        location: Option<Location>,
    ) {
        // JTS:     // if this node is already a boundary node, don't change it
        // JTS:     if (isBoundaryNode(argIndex, coord)) return;
        // if this node is already a boundary node, don't change it
        if self.is_boundary_node(arg_index, coord) {
            return;
        }

        // JTS:     if (loc == Location.BOUNDARY && useBoundaryDeterminationRule)
        // JTS:         insertBoundaryPoint(argIndex, coord);
        // JTS:     else
        // JTS:       insertPoint(argIndex, coord, loc);
        if location == Some(Location::Boundary) && self.use_boundary_determination_rule {
            self.insert_boundary_point(arg_index, coord)
        } else {
            // This assert might be overzealous (in which case we'll crash on the next line)
            // but I *think* location will never be None here (JTS doesn't really separate
            // Location.NONE but doesn't ever seem to assign it except in initializers)
            debug_assert!(location.is_some());
            self.insert_point(arg_index, coord, location.unwrap())
        }
        // JTS:   }
    }
    // JTS:
    // JTS:   // MD - experimental for now
    // JTS:   /**
    // JTS:    * Determines the {@link Location} of the given {@link Coordinate}
    // JTS:    * in this geometry.
    // JTS:    *
    // JTS:    * @param pt the point to test
    // JTS:    * @return the location of the point in the geometry
    // JTS:    */
    // JTS:   public int locate(Coordinate pt)
    // JTS:   {
    // JTS:   	if (parentGeom instanceof Polygonal && parentGeom.getNumGeometries() > 50) {
    // JTS:   		// lazily init point locator
    // JTS:   		if (areaPtLocator == null) {
    // JTS:   			areaPtLocator = new IndexedPointInAreaLocator(parentGeom);
    // JTS:   		}
    // JTS:   		return areaPtLocator.locate(pt);
    // JTS:   	}
    // JTS:   	return ptLocator.locate(pt, parentGeom);
    // JTS:   }
    // JTS: }
}
