use super::Contains;
use crate::intersects::Intersects;
use crate::{Coordinate, GeoNum, Geometry, Line, LineString, MultiPolygon, Point, Polygon};

// ┌─────────────────────────────┐
// │ Implementations for Polygon │
// └─────────────────────────────┘

impl<T> Contains<Coordinate<T>> for Polygon<T>
where
    T: GeoNum,
{
    fn contains(&self, coord: &Coordinate<T>) -> bool {
        use crate::algorithm::coordinate_position::{CoordPos, CoordinatePosition};

        self.coordinate_position(coord) == CoordPos::Inside
    }
}

impl<T> Contains<Point<T>> for Polygon<T>
where
    T: GeoNum,
{
    fn contains(&self, p: &Point<T>) -> bool {
        self.contains(&p.0)
    }
}

impl<T> Contains<Line<T>> for Polygon<T>
where
    T: GeoNum,
{
    fn contains(&self, line: &Line<T>) -> bool {
        use crate::algorithm::relate::Relate;
        Geometry::Polygon(self.clone())
            .relate(&Geometry::Line(line.clone()))
            .is_contains()
    }
}

// TODO: also check interiors
impl<T> Contains<Polygon<T>> for Polygon<T>
where
    T: GeoNum,
{
    fn contains(&self, poly: &Polygon<T>) -> bool {
        // decompose poly's exterior ring into Lines, and check each for containment
        poly.exterior().lines().all(|line| self.contains(&line))
    }
}

// TODO: ensure DE-9IM compliance
impl<T> Contains<LineString<T>> for Polygon<T>
where
    T: GeoNum,
{
    fn contains(&self, linestring: &LineString<T>) -> bool {
        // All LineString points must be inside the Polygon
        if linestring.points_iter().all(|point| self.contains(&point)) {
            // The Polygon interior is allowed to intersect with the LineString
            // but the Polygon's rings are not
            !self
                .interiors()
                .iter()
                .any(|ring| ring.intersects(linestring))
        } else {
            false
        }
    }
}

// ┌──────────────────────────────────┐
// │ Implementations for MultiPolygon │
// └──────────────────────────────────┘
// TODO: ensure DE-9IM compliance
impl<G, T> Contains<G> for MultiPolygon<T>
where
    T: CoordNum,
    Polygon<T>: Contains<G>,
{
    fn contains(&self, rhs: &G) -> bool {
        self.iter().any(|p| p.contains(rhs))
    }
}
