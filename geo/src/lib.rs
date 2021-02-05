#![doc(html_logo_url = "https://raw.githubusercontent.com/georust/meta/master/logo/logo.png")]

//! The `geo` crate provides geospatial primitive types and algorithms.
//!
//! # Types
//!
//! - **[`Coordinate`]**: A two-dimensional coordinate. All geometry types are composed of [`Coordinate`]s, though [`Coordinate`] itself is not a [`Geometry`] type.
//! - **[`Point`]**: A single point represented by one [`Coordinate`]
//! - **[`MultiPoint`]**: A collection of [`Point`]s
//! - **[`Line`]**: A line segment represented by two [`Coordinate`]s
//! - **[`LineString`]**: A series of contiguous line segments represented by two or more
//!   [`Coordinate`]s
//! - **[`MultiLineString`]**: A collection of [`LineString`]s
//! - **[`Polygon`]**: A bounded area represented by one [`LineString`] exterior ring, and zero or
//!   more [`LineString`] interior rings
//! - **[`MultiPolygon`]**: A collection of [`Polygon`]s
//! - **[`Rect`]**: An axis-aligned bounded rectangle represented by minimum and maximum
//!   [`Coordinate`]s
//! - **[`Triangle`]**: A bounded area represented by three [`Coordinate`] vertices
//! - **[`GeometryCollection`]**: A collection of [`Geometry`]s
//! - **[`Geometry`]**: An enumeration of all geometry types, excluding [`Coordinate`]
//!
//! ## Semantics
//!
//! The geospatial types provided here aim to adhere to the [OpenGIS Simple feature access][OGC-SFA]
//! standards. Thus, the types here are inter-operable with other implementations of the standards:
//! [JTS], [GEOS], etc.
//!
//! # Algorithms
//!
//! ## Area
//!
//! - **[`Area`](algorithm::area::Area)**: Calculate the planar area of a geometry
//! - **[`ChamberlainDuquetteArea`](algorithm::chamberlain_duquette_area::ChamberlainDuquetteArea)**: Calculate the geodesic area of a geometry
//!
//! ## Distance
//!
//! - **[`EuclideanDistance`](algorithm::euclidean_distance::EuclideanDistance)**: Calculate the minimum euclidean distance between geometries
//! - **[`GeodesicDistance`](algorithm::geodesic_distance::GeodesicDistance)**: Calculate the minimum geodesic distance between geometries using the algorithm presented in _Algorithms for geodesics_ by Charles Karney (2013)
//! - **[`HaversineDistance`](algorithm::haversine_distance::HaversineDistance)**: Calculate the minimum geodesic distance between geometries using the haversine formula
//! - **[`VincentyDistance`](algorithm::vincenty_distance::VincentyDistance)**: Calculate the minimum geodesic distance between geometries using Vincenty’s formula
//!
//! ## Length
//!
//! - **[`EuclideanLength`](algorithm::euclidean_length::EuclideanLength)**: Calculate the euclidean length of a geometry
//! - **[`GeodesicLength`](algorithm::geodesic_length::GeodesicLength)**: Calculate the geodesic length of a geometry using the algorithm presented in _Algorithms for geodesics_ by Charles Karney (2013)
//! - **[`HaversineLength`](algorithm::haversine_length::HaversineLength)**: Calculate the geodesic length of a geometry using the haversine formula
//! - **[`VincentyLength`](algorithm::vincenty_length::VincentyLength)**: Calculate the geodesic length of a geometry using Vincenty’s formula
//!
//! ## Simplification
//!
//! - **[`Simplify`](algorithm::simplify::Simplify)**: Simplify a geometry using the Ramer–Douglas–Peucker algorithm
//! - **[`SimplifyIdx`](algorithm::simplify::SimplifyIdx)**:
//! - **[`SimplifyVW`](algorithm::simplifyvw::SimplifyVW)**: Simplify a geometry using the Visvalingam-Whyatt algorithm
//! - **[`SimplifyVWPreserve`](algorithm::simplifyvw::SimplifyVWPreserve)**: Simplify a geometry using a topology-preserving variant of the Visvalingam-Whyatt algorithm
//! - **[`SimplifyVwIdx`](algorithm::simplifyvw::SimplifyVwIdx)**:
//!
//! ## Query
//!
//! - **[`Bearing`](algorithm::bearing::Bearing)**:
//! - **[`ClosestPoint`](algorithm::closest_point::ClosestPoint)**:
//! - **[`Contains`](algorithm::contains::Contains)**:
//! - **[`CoordinatePosition`](algorithm::coordinate_position::CoordinatePosition)**:
//! - **[`Intersects`](algorithm::intersects::Intersects)**:
//! - **[`IsConvex`](algorithm::is_convex::IsConvex)**:
//! - **[`LineInterpolatePoint`](algorithm::line_interpolate_point::LineInterpolatePoint)**:
//! - **[`LineLocatePoint`](algorithm::line_locate_point::LineLocatePoint)**:
//!
//! ## Similarity
//!
//! - **[`FrechetDistance`](algorithm::frechet_distance::FrechetDistance)**: Calculate the similarity between [`LineString`]s using the Fréchet distance
//!
//! ## Winding
//!
//! - **[`Orient`](algorithm::orient::Orient)**: Orients the exterior and interior rings of a [`Polygon`] according to convention
//! - **[`Winding`](algorithm::winding::Winding)**: Calculate and manipulate the winding order of a [`LineString`]
//!
//! ## Map
//!
//! - **[`MapCoords`](algorithm::map_coords::MapCoords)**:
//! - **[`MapCoordsInplace`](algorithm::map_coords::MapCoordsInplace)**:
//! - **[`TryMapCoords`](algorithm::map_coords::TryMapCoords)**:
//!
//! ## Boundary
//!
//! - **[`BoundingRect`](algorithm::bounding_rect::BoundingRect)**:
//! - **[`ConcaveHull`](algorithm::concave_hull::ConcaveHull)**:
//! - **[`ConvexHull`](algorithm::convex_hull::ConvexHull)**:
//! - **[`Extremes`](algorithm::extremes::Extremes)**:
//!
//! ## Affine transformations
//!
//! - **[`Rotate`](algorithm::rotate::Rotate)**
//! - **[`RotatePoint`](algorithm::rotate::RotatePoint)**:
//! - **[`Translate`](algorithm::translate::Translate)**:
//!
//! ## Unsorted
//!
//! - **[`Centroid`](algorithm::centroid::Centroid)**:
//! - **[`CoordsIter`](algorithm::coords_iter::CoordsIter)**:
//! - **[`HasDimensions`](algorithm::dimensions::HasDimensions)**:
//! - **[`HaversineDestination`](algorithm::haversine_destination::HaversineDestination)**:
//! - **[`HaversineIntermediate`](algorithm::haversine_intermediate::HaversineIntermediate)**:
//! - **[`Proj`](algorithm::proj::Proj)**:
//!
//! # Features
//!
//! The following optional [Cargo features] are available:
//!
//! - `proj-network`: Enables [network grid] support for the [`proj` crate]. After enabling this feature, [further configuration][proj crate file download] is required to use the network grid
//! - `use-proj`: Enables coordinate conversion and transformation of `Point` geometries using the [`proj` crate]
//! - `use-serde`: Allows geometry types to be serialized and deserialized with [Serde]
//!
//! [`proj` crate]: https://github.com/georust/proj
//! [Cargo features]: https://doc.rust-lang.org/cargo/reference/features.html
//! [GEOS]: https://trac.osgeo.org/geos
//! [JTS]: https://github.com/locationtech/jts
//! [network grid]: https://proj.org/usage/network.html
//! [OGC-SFA]: https://www.ogc.org/standards/sfa
//! [proj crate file download]: https://docs.rs/proj/*/proj/#grid-file-download
//! [Serde]: https://serde.rs/
//!
//! -------------------
//!
//!
//!
//!
//!
//!
//! The `geo` crate provides geospatial primitive types such as `Coordinate`, `Point`, `LineString`, and `Polygon` as
//! well as their `Multi–` equivalents, and provides algorithms and operations such as:
//!   - Area and centroid calculation
//!   - Simplification and convex hull operations
//!   - Distance measurement
//!   - Intersection checks
//!   - Affine transforms such as rotation and translation.
//!
//! The primitive types also provide the basis for other functionality in the `Geo` ecosystem, including:
//!   - Serialization to and from [GeoJSON](https://docs.rs/geojson) and [WKT](https://docs.rs/wkt)
//!   - [Coordinate transformation and projection](https://docs.rs/proj)
//!   - [Geocoding](https://docs.rs/geocoding)
//!   - [Working with GPS data](https://docs.rs/gpx)
//!
//! …allowing these crates to interoperate; GeoJSON can readily be read from a file, deserialised, transformed
//! to a local datum, modified, transformed back to `WGS84`, and serialised back to GeoJSON.
//!
//! Operations available for primitive types can be found in the `algorithm` module, along with
//! comprehensive usage examples.
//!
//! While `Geo` is primarily intended to operate on **planar** geometries, some other useful algorithms are
//! provided: Haversine, Frechet, and Vincenty distances, as well as Chamberlain-Duquette area.
//!
//!
//! ## GeoJSON
//! If you wish to read or write `GeoJSON`, use the [`geojson`](https://docs.rs/geojson) crate, with the `geo-types` feature activated.
//! This provides fallible conversions **to** `geo-types` primitives such as `Point` and `Polygon` from `geojson` `Value`
//! structs using the standard [`TryFrom`](https://doc.rust-lang.org/stable/std/convert/trait.TryFrom.html)
//! and [`TryInto`](https://doc.rust-lang.org/stable/std/convert/trait.TryInto.html) traits,
//! and conversion **from** `geo-types` primitives to `geojson`
//! `Value` structs using the [`From`](https://doc.rust-lang.org/stable/std/convert/trait.TryFrom.html) trait.

extern crate geo_types;
extern crate num_traits;
#[cfg(feature = "use-serde")]
#[macro_use]
extern crate serde;
#[cfg(feature = "use-proj")]
extern crate proj;
extern crate rstar;

pub use crate::algorithm::*;
#[allow(deprecated)]
pub use crate::traits::ToGeo;
pub use crate::types::*;

pub use geo_types::{
    line_string, point, polygon, CoordFloat, CoordNum, Coordinate, Geometry, GeometryCollection,
    Line, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon, Rect, Triangle,
};

/// This module includes all the functions of geometric calculations
pub mod algorithm;
mod traits;
mod types;
mod utils;

#[cfg(test)]
#[macro_use]
extern crate approx;

/// Mean radius of Earth in meters
/// This is the value recommended by the IUGG:
/// Moritz, H. (2000). Geodetic Reference System 1980. Journal of Geodesy, 74(1), 128–133. doi:10.1007/s001900050278
/// "Derived Geometric Constants: mean radius" (p133)
/// https://link.springer.com/article/10.1007%2Fs001900050278
/// https://sci-hub.se/https://doi.org/10.1007/s001900050278
/// https://en.wikipedia.org/wiki/Earth_radius#Mean_radius
const MEAN_EARTH_RADIUS: f64 = 6371008.8;

// Radius of Earth at the equator in meters (derived from the WGS-84 ellipsoid)
const EQUATORIAL_EARTH_RADIUS: f64 = 6_378_137.0;

// Radius of Earth at the poles in meters (derived from the WGS-84 ellipsoid)
const POLAR_EARTH_RADIUS: f64 = 6_356_752.314_245;

// Flattening of the WGS-84 ellipsoid - https://en.wikipedia.org/wiki/Flattening
const EARTH_FLATTENING: f64 =
    (EQUATORIAL_EARTH_RADIUS - POLAR_EARTH_RADIUS) / EQUATORIAL_EARTH_RADIUS;

/// A prelude which re-exports the traits for manipulating objects in this
/// crate. Typically imported with `use geo::prelude::*`.
pub mod prelude {
    pub use crate::algorithm::area::Area;
    pub use crate::algorithm::bearing::Bearing;
    pub use crate::algorithm::bounding_rect::BoundingRect;
    pub use crate::algorithm::centroid::Centroid;
    pub use crate::algorithm::chamberlain_duquette_area::ChamberlainDuquetteArea;
    pub use crate::algorithm::closest_point::ClosestPoint;
    pub use crate::algorithm::contains::Contains;
    pub use crate::algorithm::convex_hull::ConvexHull;
    pub use crate::algorithm::dimensions::HasDimensions;
    pub use crate::algorithm::euclidean_distance::EuclideanDistance;
    pub use crate::algorithm::euclidean_length::EuclideanLength;
    pub use crate::algorithm::extremes::Extremes;
    pub use crate::algorithm::frechet_distance::FrechetDistance;
    pub use crate::algorithm::geodesic_distance::GeodesicDistance;
    pub use crate::algorithm::geodesic_intermediate::GeodesicIntermediate;
    pub use crate::algorithm::geodesic_length::GeodesicLength;
    pub use crate::algorithm::haversine_destination::HaversineDestination;
    pub use crate::algorithm::haversine_distance::HaversineDistance;
    pub use crate::algorithm::haversine_intermediate::HaversineIntermediate;
    pub use crate::algorithm::haversine_length::HaversineLength;
    pub use crate::algorithm::intersects::Intersects;
    pub use crate::algorithm::is_convex::IsConvex;
    pub use crate::algorithm::map_coords::MapCoords;
    pub use crate::algorithm::orient::Orient;
    #[cfg(feature = "use-proj")]
    pub use crate::algorithm::proj::Proj;
    pub use crate::algorithm::rotate::{Rotate, RotatePoint};
    pub use crate::algorithm::simplify::Simplify;
    pub use crate::algorithm::simplifyvw::SimplifyVW;
    pub use crate::algorithm::translate::Translate;
    pub use crate::algorithm::vincenty_distance::VincentyDistance;
    pub use crate::algorithm::vincenty_length::VincentyLength;
}

/// A common numeric trait used for geo algorithms.
///
/// Different numeric types have different tradeoffs. `geo` strives to utilize generics to allow
/// users to choose their numeric types. If you are writing a function which you'd like to be
/// generic over all the numeric types supported by geo, you probably want to constraint
/// your function input to `GeoFloat`. For methods which work for integers, and not just floating
/// point, see [`GeoNum`].
///
/// # Examples
///
/// ```
/// use geo::{GeoFloat, MultiPolygon, Polygon, Point};
///
/// // An admittedly silly method implementation, but the signature shows how to use the GeoFloat trait
/// fn farthest_from<'a, T: GeoFloat>(point: &Point<T>, polygons: &'a MultiPolygon<T>) -> Option<&'a Polygon<T>> {
///     polygons.iter().fold(None, |accum, next| {
///         match accum {
///             None => Some(next),
///             Some(farthest) => {
///                 use geo::algorithm::{euclidean_distance::EuclideanDistance};
///                 if next.euclidean_distance(point) > farthest.euclidean_distance(point) {
///                     Some(next)
///                 } else {
///                     Some(farthest)
///                 }
///             }
///         }
///     })
/// }
/// ```
pub trait GeoFloat: num_traits::Float + GeoNum {}
impl<T> GeoFloat for T where T: num_traits::Float + GeoNum {}

pub trait GeoNum: CoordNum + algorithm::kernels::HasKernel {}
impl<T> GeoNum for T where T: CoordNum + algorithm::kernels::HasKernel {}
