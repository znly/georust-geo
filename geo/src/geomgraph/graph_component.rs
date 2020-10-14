use super::Label;
// JTS: import org.locationtech.jts.geom.Coordinate;
// JTS: import org.locationtech.jts.geom.IntersectionMatrix;
// JTS: import org.locationtech.jts.util.Assert;
// JTS:
// JTS: /**
// JTS:  * A GraphComponent is the parent class for the objects'
// JTS:  * that form a graph.  Each GraphComponent can carry a
// JTS:  * Label.
// JTS:  * @version 1.7
// JTS:  */
// JTS: abstract public class GraphComponent {
/// A GraphComponent is the parent class for the objects' that form a graph.
/// Each GraphComponent can carry a Label.
pub trait GraphComponent {
    // JTS:
    // JTS:   protected Label label;
    // JTS:   /**
    // JTS:    * isInResult indicates if this component has already been included in the result
    // JTS:    */
    // JTS:   private boolean isInResult = false;
    // JTS:   private boolean isCovered = false;
    // JTS:   private boolean isCoveredSet = false;
    // JTS:   private boolean isVisited = false;
    // JTS:
    // JTS:   public GraphComponent() {
    // JTS:   }
    // JTS:
    // JTS:   public GraphComponent(Label label) {
    // JTS:     this.label = label;
    // JTS:   }
    // JTS:
    // JTS:   public Label getLabel() { return label; }
    // JTS:   public void setLabel(Label label) { this.label = label; }
    fn get_label(&self) -> Option<&Label>;
    fn get_label_mut(&mut self) -> Option<&mut Label>;
    fn set_label(&mut self, new_value: Label);

    // JTS:   public void setInResult(boolean isInResult) { this.isInResult = isInResult; }
    // JTS:   public boolean isInResult() { return isInResult; }
    // JTS:   public void setCovered(boolean isCovered)
    // JTS:   {
    // JTS:     this.isCovered = isCovered;
    // JTS:     this.isCoveredSet = true;
    // JTS:   }
    // JTS:   public boolean isCovered()    { return isCovered; }
    // JTS:   public boolean isCoveredSet() { return isCoveredSet; }
    // JTS:   public boolean isVisited() { return isVisited; }
    // JTS:   public void setVisited(boolean isVisited) { this.isVisited = isVisited; }
    // JTS:   /**
    // JTS:    * @return a coordinate in this component (or null, if there are none)
    // JTS:    */
    // JTS:   abstract public Coordinate getCoordinate();
    // JTS:   /**
    // JTS:    * compute the contribution to an IM for this component
    // JTS:    */
    // JTS:   abstract protected void computeIM(IntersectionMatrix im);
    // JTS:   /**
    // JTS:    * An isolated component is one that does not intersect or touch any other
    // JTS:    * component.  This is the case if the label has valid locations for
    // JTS:    * only a single Geometry.
    // JTS:    *
    // JTS:    * @return true if this component is isolated
    // JTS:    */
    // JTS:   abstract public boolean isIsolated();
    // JTS:   /**
    // JTS:    * Update the IM with the contribution for this component.
    // JTS:    * A component only contributes if it has a labelling for both parent geometries
    // JTS:    */
    // JTS:   public void updateIM(IntersectionMatrix im)
    // JTS:   {
    // JTS:     Assert.isTrue(label.getGeometryCount() >= 2, "found partial label");
    // JTS:     computeIM(im);
    // JTS:   }
    // JTS:
    // JTS: }
}
