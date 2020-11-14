use crate::{
    bounding_rect::BoundingRect,
    coordinate_position::{CoordPos, CoordinatePosition},
    dimensions::HasDimensions,
    geomgraph::{
        algorithm::boundary_node_rule::{BoundaryNodeRule, Mod2BoundaryNodeRule},
        Float, Location,
    },
    intersects::Intersects,
    Coordinate, Geometry, Line, LineString, Point, Polygon,
};
// JTS: /**
// JTS:  * Computes the topological ({@link Location})
// JTS:  * of a single point to a {@link Geometry}.
// JTS:  * A {@link BoundaryNodeRule} may be specified
// JTS:  * to control the evaluation of whether the point lies on the boundary or not
// JTS:  * The default rule is to use the the <i>SFS Boundary Determination Rule</i>
// JTS:  * <p>
// JTS:  * Notes:
// JTS:  * <ul>
// JTS:  * <li>{@link LinearRing}s do not enclose any area - points inside the ring are still in the EXTERIOR of the ring.
// JTS:  * </ul>
// JTS:  * Instances of this class are not reentrant.
// JTS:  *
// JTS:  * @version 1.7
// JTS:  */
// JTS: public class PointLocator
// JTS: {
pub(crate) struct PointLocator<F>
where
    F: Float,
{
    is_in: bool,
    num_boundaries: usize,
    _marker: std::marker::PhantomData<F>,
}

impl<F> PointLocator<F>
where
    F: Float,
{
    // JTS:   // default is to use OGC SFS rule
    // JTS:   private BoundaryNodeRule boundaryRule =
    // JTS:   	//BoundaryNodeRule.ENDPOINT_BOUNDARY_RULE;
    // JTS:   	BoundaryNodeRule.OGC_SFS_BOUNDARY_RULE;
    // JTS:
    // JTS:   private boolean isIn;         // true if the point lies in or on any Geometry element
    // JTS:   private int numBoundaries;    // the number of sub-elements whose boundaries the point lies in
    // JTS:
    // JTS:   public PointLocator() {
    // JTS:   }
    pub fn new() -> Self {
        PointLocator {
            is_in: false,
            num_boundaries: 0,
            _marker: std::marker::PhantomData,
        }
    }

    // JTS:   public PointLocator(BoundaryNodeRule boundaryRule)
    // JTS:   {
    // JTS:     if (boundaryRule == null)
    // JTS:       throw new IllegalArgumentException("Rule must be non-null");
    // JTS:     this.boundaryRule = boundaryRule;
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Convenience method to test a point for intersection with
    // JTS:    * a Geometry
    // JTS:    * @param p the coordinate to test
    // JTS:    * @param geom the Geometry to test
    // JTS:    * @return <code>true</code> if the point is in the interior or boundary of the Geometry
    // JTS:    */
    // JTS:   public boolean intersects(Coordinate p, Geometry geom)
    // JTS:   {
    // JTS:     return locate(p, geom) != Location.EXTERIOR;
    // JTS:   }
    // JTS:
    // JTS:   /**
    // JTS:    * Computes the topological relationship ({@link Location}) of a single point
    // JTS:    * to a Geometry.
    // JTS:    * It handles both single-element
    // JTS:    * and multi-element Geometries.
    // JTS:    * The algorithm for multi-part Geometries
    // JTS:    * takes into account the SFS Boundary Determination Rule.
    // JTS:    *
    // JTS:    * @return the {@link Location} of the point relative to the input Geometry
    // JTS:    */
    // JTS:   public int locate(Coordinate p, Geometry geom)
    // JTS:   {
    pub fn locate(&self, coordinate: &Coordinate<F>, geometry: &Geometry<F>) -> Location {
        match geometry.coordinate_position(coordinate) {
            CoordPos::OnBoundary => Location::Boundary,
            CoordPos::Inside => Location::Interior,
            CoordPos::Outside => Location::Exterior,
        }
    }

    //     // CLEANUP: remove this class - all the functionality lives in the new coordinate_position module
    //     // JTS:     if (geom.isEmpty()) return Location.EXTERIOR;
    //     if geometry.is_empty() {
    //         return Location::Exterior;
    //     }
    //
    //     match geometry {
    //         Geometry::LineString(line_string) => {
    //             // JTS:     if (geom instanceof LineString) {
    //             // JTS:       return locateOnLineString(p, (LineString) geom);
    //             // JTS:     }
    //             self.locate_on_line_string(coordinate, line_string)
    //         }
    //         Geometry::Polygon(polygon) => {
    //             // JTS:     else if (geom instanceof Polygon) {
    //             // JTS:       return locateInPolygon(p, (Polygon) geom);
    //             // JTS:     }
    //             self.locate_in_polygon(coordinate, polygon)
    //         }
    //         _ => {
    //             // JTS:     isIn = false;
    //             // JTS:     numBoundaries = 0;
    //             // JTS:     computeLocation(p, geom);
    //             // JTS:     if (boundaryRule.isInBoundary(numBoundaries))
    //             // JTS:       return Location.BOUNDARY;
    //             // JTS:     if (numBoundaries > 0 || isIn)
    //             // JTS:       return Location.INTERIOR;
    //             // JTS:
    //             // JTS:     return Location.EXTERIOR;
    //             // JTS:   }
    //             self.is_in = false;
    //             self.num_boundaries = 0;
    //             self.compute_location(coordinate, geometry);
    //
    //             if Mod2BoundaryNodeRule.is_in_boundary(self.num_boundaries) {
    //                 Location::Boundary
    //             } else if self.num_boundaries > 0 || self.is_in {
    //                 Location::Interior
    //             } else {
    //                 Location::Exterior
    //             }
    //         }
    //     }
    // }
    //
    // // JTS:   private void computeLocation(Coordinate p, Geometry geom)
    // // JTS:   {
    // fn compute_location(&mut self, coordinate: &Coordinate<F>, geometry: &Geometry<F>) {
    //     match geometry {
    //         // JTS:     if (geom instanceof Point) {
    //         // JTS:       updateLocationInfo(locateOnPoint(p, (Point) geom));
    //         // JTS:     }
    //         Geometry::Point(point) => {
    //             self.update_location_info(self.locate_on_point(coordinate, point))
    //         }
    //         // JTS:     if (geom instanceof LineString) {
    //         // JTS:       updateLocationInfo(locateOnLineString(p, (LineString) geom));
    //         // JTS:     }
    //         Geometry::LineString(line_string) => {
    //             self.update_location_info(self.locate_on_line_string(coordinate, line_string))
    //         }
    //         // JTS:     else if (geom instanceof Polygon) {
    //         // JTS:       updateLocationInfo(locateInPolygon(p, (Polygon) geom));
    //         // JTS:     }
    //         Geometry::Polygon(polygon) => {
    //             self.update_location_info(self.locate_in_polygon(coordinate, polygon))
    //         }
    //         // JTS:     else if (geom instanceof MultiLineString) {
    //         // JTS:       MultiLineString ml = (MultiLineString) geom;
    //         // JTS:       for (int i = 0; i < ml.getNumGeometries(); i++) {
    //         // JTS:         LineString l = (LineString) ml.getGeometryN(i);
    //         // JTS:         updateLocationInfo(locateOnLineString(p, l));
    //         // JTS:       }
    //         // JTS:     }
    //         Geometry::MultiLineString(multi_line_string) => {
    //             for line_string in &multi_line_string.0 {
    //                 self.update_location_info(self.locate_on_line_string(coordinate, line_string));
    //             }
    //         }
    //         // JTS:     else if (geom instanceof MultiPolygon) {
    //         // JTS:       MultiPolygon mpoly = (MultiPolygon) geom;
    //         // JTS:       for (int i = 0; i < mpoly.getNumGeometries(); i++) {
    //         // JTS:         Polygon poly = (Polygon) mpoly.getGeometryN(i);
    //         // JTS:         updateLocationInfo(locateInPolygon(p, poly));
    //         // JTS:       }
    //         // JTS:     }
    //         Geometry::MultiPolygon(multi_polygon) => {
    //             for polygon in &multi_polygon.0 {
    //                 self.update_location_info(self.locate_in_polygon(coordinate, polygon))
    //             }
    //         }
    //         // JTS:     else if (geom instanceof GeometryCollection) {
    //         // JTS:       Iterator geomi = new GeometryCollectionIterator((GeometryCollection) geom);
    //         // JTS:       while (geomi.hasNext()) {
    //         // JTS:         Geometry g2 = (Geometry) geomi.next();
    //         // JTS:         if (g2 != geom)
    //         // JTS:           computeLocation(p, g2);
    //         // JTS:       }
    //         // JTS:     }
    //         // JTS:   }
    //         Geometry::GeometryCollection(geometry_collection) => {
    //             for geometry in geometry_collection {
    //                 self.compute_location(coordinate, geometry)
    //             }
    //         }
    //         Geometry::Line(line) => {
    //             self.update_location_info(self.locate_on_line(coordinate, line))
    //         }
    //         Geometry::MultiPoint(multi_point) => {
    //             for point in &multi_point.0 {
    //                 self.update_location_info(self.locate_on_point(coordinate, point))
    //             }
    //         }
    //         Geometry::Rect(rect) => {
    //             // PERF: avoid clone, implement Rect specific handling
    //             self.update_location_info(
    //                 self.locate_in_polygon(coordinate, &Polygon::from(rect.clone())),
    //             );
    //         }
    //         Geometry::Triangle(triangle) => {
    //             // PERF: avoid clone, implement Triangle specific handling
    //             self.update_location_info(
    //                 self.locate_in_polygon(coordinate, &Polygon::from(triangle.clone())),
    //             );
    //         }
    //     }
    // }
    // // JTS:
    // // JTS:   private void updateLocationInfo(int loc)
    // // JTS:   {
    // // JTS:     if (loc == Location.INTERIOR) isIn = true;
    // // JTS:     if (loc == Location.BOUNDARY) numBoundaries++;
    // // JTS:   }
    // fn update_location_info(&mut self, location: Location) {
    //     match location {
    //         Location::Interior => {
    //             self.is_in = true;
    //         }
    //         Location::Boundary => {
    //             self.num_boundaries += 1;
    //         }
    //         Location::Exterior => {}
    //     }
    // }
    //
    // // JTS:   private int locateOnPoint(Coordinate p, Point pt)
    // // JTS:   {
    // fn locate_on_point(&self, coordinate: &Coordinate<F>, point: &Point<F>) -> Location {
    //     // JTS:   	// no point in doing envelope test, since equality test is just as fast
    //     // JTS:
    //     // JTS:     Coordinate ptCoord = pt.getCoordinate();
    //     // JTS:     if (ptCoord.equals2D(p))
    //     // JTS:       return Location.INTERIOR;
    //     // JTS:     return Location.EXTERIOR;
    //     // JTS:   }
    //     if &Point::from(*coordinate) == point {
    //         Location::Interior
    //     } else {
    //         Location::Exterior
    //     }
    // }
    //
    // // JTS:   private int locateOnLineString(Coordinate p, LineString l)
    // // JTS:   {
    // fn locate_on_line_string(
    //     &self,
    //     coordinate: &Coordinate<F>,
    //     line_string: &LineString<F>,
    // ) -> Location {
    //     if line_string.0.len() < 2 {
    //         debug_assert!(false, "invalid line string with less than 2 points");
    //         return Location::Exterior;
    //     }
    //
    //     // JTS:     // bounding-box check
    //     // JTS:     if (! l.getEnvelopeInternal().intersects(p)) return Location.EXTERIOR;
    //     if !line_string.bounding_rect().unwrap().intersects(coordinate) {
    //         return Location::Exterior;
    //     }
    //
    //     // JTS:     CoordinateSequence seq = l.getCoordinateSequence();
    //     // JTS:     if (! l.isClosed()) {
    //     // JTS:           if (p.equals(seq.getCoordinate(0))
    //     // JTS:           || p.equals(seq.getCoordinate(seq.size() - 1)) ) {
    //     // JTS:         return Location.BOUNDARY;
    //     // JTS:       }
    //     // JTS:     }
    //     if !line_string.is_closed() {
    //         if coordinate == line_string.0.first().unwrap()
    //             || coordinate == line_string.0.last().unwrap()
    //         {
    //             return Location::Boundary;
    //         }
    //     }
    //
    //     // JTS:     if (PointLocation.isOnLine(p, seq)) {
    //     // JTS:       return Location.INTERIOR;
    //     // JTS:     }
    //     // JTS:     return Location.EXTERIOR;
    //     // JTS:   }
    //     if line_string.intersects(coordinate) {
    //         // We've already checked for "Boundary" condition, so if there's an intersection at
    //         // this point, it must be on `line_string`'s interior
    //         Location::Interior
    //     } else {
    //         Location::Exterior
    //     }
    // }
    //
    // fn locate_on_line(&self, coordinate: &Coordinate<F>, line: &Line<F>) -> Location {
    //     // degenerate line is a point
    //     if line.start == line.end {
    //         if &line.start == coordinate {
    //             return Location::Interior;
    //         } else {
    //             return Location::Exterior;
    //         }
    //     }
    //
    //     if coordinate == &line.start || coordinate == &line.end {
    //         return Location::Boundary;
    //     } else if line.intersects(coordinate) {
    //         Location::Interior
    //     } else {
    //         Location::Exterior
    //     }
    // }
    //
    // // JTS:   private int locateInPolygonRing(Coordinate p, LinearRing ring)
    // // JTS:   {
    // // JTS:   	// bounding-box check
    // // JTS:   	if (! ring.getEnvelopeInternal().intersects(p)) return Location.EXTERIOR;
    // // JTS:
    // // JTS:   	return PointLocation.locateInRing(p, ring.getCoordinates());
    // // JTS:   }
    // fn locate_in_polygon_ring(&self, coordinate: &Coordinate<F>, ring: &LineString<F>) -> Location {
    //     // CLEANUP: remove unwrap?
    //     if !ring.bounding_rect().unwrap().intersects(coordinate) {
    //         return Location::Exterior;
    //     }
    //
    //     use crate::utils::coord_pos_relative_to_ring;
    //     coord_pos_relative_to_ring(*coordinate, ring).into()
    // }
    //
    // // JTS:   private int locateInPolygon(Coordinate p, Polygon poly)
    // // JTS:   {
    // fn locate_in_polygon(&self, coordinate: &Coordinate<F>, polygon: &Polygon<F>) -> Location {
    //     // JTS:     if (poly.isEmpty()) return Location.EXTERIOR;
    //     if polygon.is_empty() {
    //         return Location::Exterior;
    //     }
    //
    //     // JTS:     LinearRing shell = poly.getExteriorRing();
    //     // JTS:     int shellLoc = locateInPolygonRing(p, shell);
    //     let shell_location = self.locate_in_polygon_ring(coordinate, polygon.exterior());
    //     // JTS:     if (shellLoc == Location.EXTERIOR) return Location.EXTERIOR;
    //     // JTS:     if (shellLoc == Location.BOUNDARY) return Location.BOUNDARY;
    //     if shell_location == Location::Exterior {
    //         return Location::Exterior;
    //     }
    //     if shell_location == Location::Boundary {
    //         return Location::Boundary;
    //     }
    //
    //     debug_assert!(shell_location == Location::Interior);
    //     // JTS:     // now test if the point lies in or on the holes
    //     // JTS:     for (int i = 0; i < poly.getNumInteriorRing(); i++) {
    //     // JTS:       LinearRing hole = poly.getInteriorRingN(i);
    //     // JTS:       int holeLoc = locateInPolygonRing(p, hole);
    //     // JTS:       if (holeLoc == Location.INTERIOR) return Location.EXTERIOR;
    //     // JTS:       if (holeLoc == Location.BOUNDARY) return Location.BOUNDARY;
    //     // JTS:     }
    //     for hole in polygon.interiors() {
    //         let hole_location = self.locate_in_polygon_ring(coordinate, hole);
    //         if hole_location == Location::Interior {
    //             return Location::Exterior;
    //         }
    //         if hole_location == Location::Boundary {
    //             return Location::Boundary;
    //         }
    //     }
    //
    //     // JTS:     return Location.INTERIOR;
    //     // JTS:   }
    //     Location::Interior
    // }
    // JTS: }
}
