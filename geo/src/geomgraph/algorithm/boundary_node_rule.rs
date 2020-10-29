// JTS: package org.locationtech.jts.algorithm;
// JTS:
// JTS: import org.locationtech.jts.geom.LineString;
// JTS: import org.locationtech.jts.geom.Lineal;
// JTS: import org.locationtech.jts.geom.LinearRing;
// JTS: import org.locationtech.jts.geom.MultiLineString;
// JTS: import org.locationtech.jts.operation.BoundaryOp;
// JTS: import org.locationtech.jts.operation.IsSimpleOp;
// JTS: import org.locationtech.jts.operation.relate.RelateOp;
// JTS:
// JTS: /**
// JTS:  * An interface for rules which determine whether node points
// JTS:  * which are in boundaries of {@link Lineal} geometry components
// JTS:  * are in the boundary of the parent geometry collection.
// JTS:  * The SFS specifies a single kind of boundary node rule,
// JTS:  * the {@link Mod2BoundaryNodeRule} rule.
// JTS:  * However, other kinds of Boundary Node Rules are appropriate
// JTS:  * in specific situations (for instance, linear network topology
// JTS:  * usually follows the {@link EndPointBoundaryNodeRule}.)
// JTS:  * Some JTS operations
// JTS:  * (such as {@link RelateOp}, {@link BoundaryOp} and {@link IsSimpleOp})
// JTS:  * allow the BoundaryNodeRule to be specified,
// JTS:  * and respect the supplied rule when computing the results of the operation.
// JTS:  * <p>
// JTS:  * An example use case for a non-SFS-standard Boundary Node Rule is
// JTS:  * that of checking that a set of {@link LineString}s have
// JTS:  * valid linear network topology, when turn-arounds are represented
// JTS:  * as closed rings.  In this situation, the entry road to the
// JTS:  * turn-around is only valid when it touches the turn-around ring
// JTS:  * at the single (common) endpoint.  This is equivalent
// JTS:  * to requiring the set of <tt>LineString</tt>s to be
// JTS:  * <b>simple</b> under the {@link EndPointBoundaryNodeRule}.
// JTS:  * The SFS-standard {@link Mod2BoundaryNodeRule} is not
// JTS:  * sufficient to perform this test, since it
// JTS:  * states that closed rings have <b>no</b> boundary points.
// JTS:  * <p>
// JTS:  * This interface and its subclasses follow the <tt>Strategy</tt> design pattern.
// JTS:  *
// JTS:  * @author Martin Davis
// JTS:  * @version 1.7
// JTS:  *
// JTS:  * @see RelateOp
// JTS:  * @see BoundaryOp
// JTS:  * @see IsSimpleOp
// JTS:  * @see PointLocator
// JTS:  */
// JTS: public interface BoundaryNodeRule
// JTS: {
// JTS:
// JTS: 	/**
// JTS: 	 * Tests whether a point that lies in <tt>boundaryCount</tt>
// JTS: 	 * geometry component boundaries is considered to form part of the boundary
// JTS: 	 * of the parent geometry.
// JTS: 	 *
// JTS: 	 * @param boundaryCount the number of component boundaries that this point occurs in
// JTS: 	 * @return true if points in this number of boundaries lie in the parent boundary
// JTS: 	 */
// JTS:   boolean isInBoundary(int boundaryCount);
pub trait BoundaryNodeRule {
    fn is_in_boundary(&self, boundary_count: usize) -> bool;
}

// JTS:
// JTS:   /**
// JTS:    * The Mod-2 Boundary Node Rule (which is the rule specified in the OGC SFS).
// JTS:    * @see Mod2BoundaryNodeRule
// JTS:    */
// JTS:   public static final BoundaryNodeRule MOD2_BOUNDARY_RULE = new Mod2BoundaryNodeRule();
// JTS:
// JTS:   /**
// JTS:    * The Endpoint Boundary Node Rule.
// JTS:    * @see EndPointBoundaryNodeRule
// JTS:    */
// JTS:   public static final BoundaryNodeRule ENDPOINT_BOUNDARY_RULE = new EndPointBoundaryNodeRule();
// JTS:
// JTS:   /**
// JTS:    * The MultiValent Endpoint Boundary Node Rule.
// JTS:    * @see MultiValentEndPointBoundaryNodeRule
// JTS:    */
// JTS:   public static final BoundaryNodeRule MULTIVALENT_ENDPOINT_BOUNDARY_RULE = new MultiValentEndPointBoundaryNodeRule();
// JTS:
// JTS:   /**
// JTS:    * The Monovalent Endpoint Boundary Node Rule.
// JTS:    * @see MonoValentEndPointBoundaryNodeRule
// JTS:    */
// JTS:   public static final BoundaryNodeRule MONOVALENT_ENDPOINT_BOUNDARY_RULE = new MonoValentEndPointBoundaryNodeRule();
// JTS:
// JTS:   /**
// JTS:    * The Boundary Node Rule specified by the OGC Simple Features Specification,
// JTS:    * which is the same as the Mod-2 rule.
// JTS:    * @see Mod2BoundaryNodeRule
// JTS:    */
// JTS:   public static final BoundaryNodeRule OGC_SFS_BOUNDARY_RULE = MOD2_BOUNDARY_RULE;
// JTS:
// JTS:   /**
// JTS:    * A {@link BoundaryNodeRule} specifies that points are in the
// JTS:    * boundary of a lineal geometry iff
// JTS:    * the point lies on the boundary of an odd number
// JTS:    * of components.
// JTS:    * Under this rule {@link LinearRing}s and closed
// JTS:    * {@link LineString}s have an empty boundary.
// JTS:    * <p>
// JTS:    * This is the rule specified by the <i>OGC SFS</i>,
// JTS:    * and is the default rule used in JTS.
// JTS:    *
// JTS:    * @author Martin Davis
// JTS:    * @version 1.7
// JTS:    */
// JTS:   public static class Mod2BoundaryNodeRule
// JTS:       implements BoundaryNodeRule
// JTS:   {
// JTS:     public boolean isInBoundary(int boundaryCount)
// JTS:     {
// JTS:       // the "Mod-2 Rule"
// JTS:       return boundaryCount % 2 == 1;
// JTS:     }
// JTS:   }
/// A [BoundaryNodeRule]() specifies that points are in the boundary of a
/// lineal geometry iff the point lies on the boundary of an odd number of
/// components.
///
/// Under this rule LinearRings and closed LineStrings have an empty boundary.
///
/// This is the rule specified by the _OGC SFS_, and is the default rule
pub(crate) struct Mod2BoundaryNodeRule;
impl BoundaryNodeRule for Mod2BoundaryNodeRule {
    fn is_in_boundary(&self, boundary_count: usize) -> bool {
        // the "Mod-2 Rule"
        return boundary_count % 2 == 1;
    }
}

// JTS:   /**
// JTS:    * A {@link BoundaryNodeRule} which specifies that any points which are endpoints
// JTS:    * of lineal components are in the boundary of the
// JTS:    * parent geometry.
// JTS:    * This corresponds to the "intuitive" topological definition
// JTS:    * of boundary.
// JTS:    * Under this rule {@link LinearRing}s have a non-empty boundary
// JTS:    * (the common endpoint of the underlying LineString).
// JTS:    * <p>
// JTS:    * This rule is useful when dealing with linear networks.
// JTS:    * For example, it can be used to check
// JTS:    * whether linear networks are correctly noded.
// JTS:    * The usual network topology constraint is that linear segments may touch only at endpoints.
// JTS:    * In the case of a segment touching a closed segment (ring) at one point,
// JTS:    * the Mod2 rule cannot distinguish between the permitted case of touching at the
// JTS:    * node point and the invalid case of touching at some other interior (non-node) point.
// JTS:    * The EndPoint rule does distinguish between these cases,
// JTS:    * so is more appropriate for use.
// JTS:    *
// JTS:    * @author Martin Davis
// JTS:    * @version 1.7
// JTS:    */
// JTS:   public static class EndPointBoundaryNodeRule
// JTS:       implements BoundaryNodeRule
// JTS:   {
// JTS:     public boolean isInBoundary(int boundaryCount)
// JTS:     {
// JTS:       return boundaryCount > 0;
// JTS:     }
// JTS:   }
// JTS:
// JTS:   /**
// JTS:    * A {@link BoundaryNodeRule} which determines that only
// JTS:    * endpoints with valency greater than 1 are on the boundary.
// JTS:    * This corresponds to the boundary of a {@link MultiLineString}
// JTS:    * being all the "attached" endpoints, but not
// JTS:    * the "unattached" ones.
// JTS:    *
// JTS:    * @author Martin Davis
// JTS:    * @version 1.7
// JTS:    */
// JTS:   public static class MultiValentEndPointBoundaryNodeRule
// JTS:       implements BoundaryNodeRule
// JTS:   {
// JTS:     public boolean isInBoundary(int boundaryCount)
// JTS:     {
// JTS:       return boundaryCount > 1;
// JTS:     }
// JTS:   }
// JTS:
// JTS:   /**
// JTS:    * A {@link BoundaryNodeRule} which determines that only
// JTS:    * endpoints with valency of exactly 1 are on the boundary.
// JTS:    * This corresponds to the boundary of a {@link MultiLineString}
// JTS:    * being all the "unattached" endpoints.
// JTS:    *
// JTS:    * @author Martin Davis
// JTS:    * @version 1.7
// JTS:    */
// JTS:   public static class MonoValentEndPointBoundaryNodeRule
// JTS:       implements BoundaryNodeRule
// JTS:   {
// JTS:     public boolean isInBoundary(int boundaryCount)
// JTS:     {
// JTS:       return boundaryCount == 1;
// JTS:     }
// JTS:   }
// JTS:
// JTS:
// JTS: }
