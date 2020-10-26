// JTS: import org.locationtech.jts.geom.Location;
use super::{Location, Position, TopologyLocation};

// JTS:
// JTS:  /**
// JTS:  * A <code>Label</code> indicates the topological relationship of a component
// JTS:  * of a topology graph to a given <code>Geometry</code>.
// JTS:  * This class supports labels for relationships to two <code>Geometry</code>s,
// JTS:  * which is sufficient for algorithms for binary operations.
// JTS:  * <P>
// JTS:  * Topology graphs support the concept of labeling nodes and edges in the graph.
// JTS:  * The label of a node or edge specifies its topological relationship to one or
// JTS:  * more geometries.  (In fact, since JTS operations have only two arguments labels
// JTS:  * are required for only two geometries).  A label for a node or edge has one or
// JTS:  * two elements, depending on whether the node or edge occurs in one or both of the
// JTS:  * input <code>Geometry</code>s.  Elements contain attributes which categorize the
// JTS:  * topological location of the node or edge relative to the parent
// JTS:  * <code>Geometry</code>; that is, whether the node or edge is in the interior,
// JTS:  * boundary or exterior of the <code>Geometry</code>.  Attributes have a value
// JTS:  * from the set <code>{Interior, Boundary, Exterior}</code>.  In a node each
// JTS:  * element has  a single attribute <code>&lt;On&gt;</code>.  For an edge each element has a
// JTS:  * triplet of attributes <code>&lt;Left, On, Right&gt;</code>.
// JTS:  * <P>
// JTS:  * It is up to the client code to associate the 0 and 1 <code>TopologyLocation</code>s
// JTS:  * with specific geometries.
// JTS:  * @version 1.7
// JTS:  *
// JTS:  */
// JTS: public class Label {
#[derive(Clone)]
pub struct Label {
    // REVIEW: better name? what does this stand for - "element location's topology"?
    elt: [TopologyLocation; 2],
}

impl Label {
    // JTS:
    // JTS:   // converts a Label to a Line label (that is, one with no side Locations)
    // JTS:   public static Label toLineLabel(Label label)
    // JTS:   {
    // JTS:     Label lineLabel = new Label(Location.NONE);
    // JTS:     for (int i = 0; i < 2; i++) {
    // JTS:       lineLabel.setLocation(i, label.getLocation(i));
    // JTS:     }
    // JTS:     return lineLabel;
    // JTS:   }
    // JTS:
    // JTS:   TopologyLocation elt[] = new TopologyLocation[2];
    // JTS:
    // JTS:   /**
    // JTS:    * Construct a Label with a single location for both Geometries.
    // JTS:    * Initialize the locations to Null
    // JTS:    */
    // JTS:   public Label(int onLoc)
    // JTS:   {
    // JTS:     elt[0] = new TopologyLocation(onLoc);
    // JTS:     elt[1] = new TopologyLocation(onLoc);
    // JTS:   }
    // JTS:   /**
    // JTS:    * Construct a Label with a single location for both Geometries.
    // JTS:    * Initialize the location for the Geometry index.
    // JTS:    */
    // JTS:   public Label(int geomIndex, int onLoc)
    // JTS:   {
    // JTS:     elt[0] = new TopologyLocation(Location.NONE);
    // JTS:     elt[1] = new TopologyLocation(Location.NONE);
    // JTS:     elt[geomIndex].setLocation(onLoc);
    // JTS:   }
    pub fn new_with_on_location(geom_index: usize, location: Option<Location>) -> Label {
        let mut label = Label {
            elt: [
                TopologyLocation::new_on(None),
                TopologyLocation::new_on(None),
            ],
        };
        if let Some(location) = location {
            label.elt[geom_index].set_on_location(location);
        }

        label
    }

    // JTS:   /**
    // JTS:    * Construct a Label with On, Left and Right locations for both Geometries.
    // JTS:    * Initialize the locations for both Geometries to the given values.
    // JTS:    */
    // JTS:   public Label(int onLoc, int leftLoc, int rightLoc)
    // JTS:   {
    // JTS:     elt[0] = new TopologyLocation(onLoc, leftLoc, rightLoc);
    // JTS:     elt[1] = new TopologyLocation(onLoc, leftLoc, rightLoc);
    // JTS:   }
    // JTS:   /**
    // JTS:    * Construct a Label with On, Left and Right locations for both Geometries.
    // JTS:    * Initialize the locations for the given Geometry index.
    // JTS:    */
    // JTS:   public Label(int geomIndex, int onLoc, int leftLoc, int rightLoc)
    // JTS:   {
    // JTS:     elt[0] = new TopologyLocation(Location.NONE, Location.NONE, Location.NONE);
    // JTS:     elt[1] = new TopologyLocation(Location.NONE, Location.NONE, Location.NONE);
    // JTS:     elt[geomIndex].setLocations(onLoc, leftLoc, rightLoc);
    // JTS:   }
    pub fn new_with_locations(
        geom_index: usize,
        on_location: Location,
        left_location: Location,
        right_location: Location,
    ) -> Label {
        let mut label = Label {
            elt: [
                TopologyLocation::new_on_left_right(None, None, None),
                TopologyLocation::new_on_left_right(None, None, None),
            ],
        };
        label.elt[geom_index].set_locations(on_location, left_location, right_location);
        label
    }
    // JTS:   /**
    // JTS:    * Construct a Label with the same values as the argument Label.
    // JTS:    */
    // JTS:   public Label(Label lbl)
    // JTS:   {
    // JTS:     elt[0] = new TopologyLocation(lbl.elt[0]);
    // JTS:     elt[1] = new TopologyLocation(lbl.elt[1]);
    // JTS:   }
    // JTS:
    // JTS:   public void flip()
    // JTS:   {
    // JTS:     elt[0].flip();
    // JTS:     elt[1].flip();
    // JTS:   }

    // JTS:   public int getLocation(int geomIndex, int posIndex) { return elt[geomIndex].get(posIndex); }
    pub fn location(&self, geom_index: usize, position: Position) -> Option<Location> {
        self.elt[geom_index].get(position)
    }

    // JTS:   public int getLocation(int geomIndex) { return elt[geomIndex].get(Position.ON); }
    pub fn on_location(&self, geom_index: usize) -> Option<Location> {
        self.elt[geom_index].get(Position::On)
    }

    // JTS:   public void setLocation(int geomIndex, int posIndex, int location)
    // JTS:   {
    // JTS:     elt[geomIndex].setLocation(posIndex, location);
    // JTS:   }
    pub fn set_location(&mut self, geom_index: usize, position: Position, location: Location) {
        self.elt[geom_index].set_location(position, location);
    }

    // JTS:   public void setLocation(int geomIndex, int location)
    // JTS:   {
    // JTS:     elt[geomIndex].setLocation(Position.ON, location);
    // JTS:   }
    pub fn set_on_location(&mut self, geom_index: usize, location: Location) {
        self.elt[geom_index].set_location(Position::On, location);
    }

    // JTS:   public void setAllLocations(int geomIndex, int location)
    // JTS:   {
    // JTS:     elt[geomIndex].setAllLocations(location);
    // JTS:   }
    pub fn set_all_locations(&mut self, geom_index: usize, location: Location) {
        self.elt[geom_index].set_all_locations(location)
    }

    // JTS:   public void setAllLocationsIfNull(int geomIndex, int location)
    // JTS:   {
    // JTS:     elt[geomIndex].setAllLocationsIfNull(location);
    // JTS:   }
    // JTS:   public void setAllLocationsIfNull(int location)
    // JTS:   {
    // JTS:     setAllLocationsIfNull(0, location);
    // JTS:     setAllLocationsIfNull(1, location);
    // JTS:   }
    // JTS:   /**
    // JTS:    * Merge this label with another one.
    // JTS:    * Merging updates any null attributes of this label with the attributes from lbl
    // JTS:    */
    // JTS:   public void merge(Label lbl)
    // JTS:   {
    // JTS:     for (int i = 0; i < 2; i++) {
    // JTS:       if (elt[i] == null && lbl.elt[i] != null) {
    // JTS:         elt[i] = new TopologyLocation(lbl.elt[i]);
    // JTS:       }
    // JTS:       else {
    // JTS:         elt[i].merge(lbl.elt[i]);
    // JTS:       }
    // JTS:     }
    // JTS:   }
    // JTS:   public int getGeometryCount()
    // JTS:   {
    // JTS:     int count = 0;
    // JTS:     if (! elt[0].isNull()) count++;
    // JTS:     if (! elt[1].isNull()) count++;
    // JTS:     return count;
    // JTS:   }
    pub fn geometry_count(&self) -> usize {
        self.elt
            .iter()
            .filter(|location| !location.is_empty())
            .count()
    }

    // JTS:   public boolean isNull(int geomIndex) { return elt[geomIndex].isNull(); }
    pub fn is_empty(&self, geom_index: usize) -> bool {
        self.elt[geom_index].is_empty()
    }

    // JTS:   public boolean isAnyNull(int geomIndex) { return elt[geomIndex].isAnyNull(); }
    // JTS:
    // JTS:   public boolean isArea()               { return elt[0].isArea() || elt[1].isArea();   }
    // JTS:   public boolean isArea(int geomIndex)
    // JTS:   {
    // JTS:   	/*  Testing
    // JTS:   	if (elt[0].getLocations().length != elt[1].getLocations().length) {
    // JTS:   		System.out.println(this);
    // JTS:   	}
    // JTS:   		*/
    // JTS:   	return elt[geomIndex].isArea();
    // JTS:   }
    // JTS:   public boolean isLine(int geomIndex)  { return elt[geomIndex].isLine();   }
    // JTS:
    // JTS:   public boolean isEqualOnSide(Label lbl, int side)
    // JTS:   {
    // JTS:     return
    // JTS:           this.elt[0].isEqualOnSide(lbl.elt[0], side)
    // JTS:       &&  this.elt[1].isEqualOnSide(lbl.elt[1], side);
    // JTS:   }
    // JTS:   public boolean allPositionsEqual(int geomIndex, int loc)
    // JTS:   {
    // JTS:     return elt[geomIndex].allPositionsEqual(loc);
    // JTS:   }
    // JTS:   /**
    // JTS:    * Converts one GeometryLocation to a Line location
    // JTS:    */
    // JTS:   public void toLine(int geomIndex)
    // JTS:   {
    // JTS:     if (elt[geomIndex].isArea())
    // JTS:       elt[geomIndex] = new TopologyLocation(elt[geomIndex].location[0]);
    // JTS:   }
    // JTS:   public String toString()
    // JTS:   {
    // JTS:     StringBuffer buf = new StringBuffer();
    // JTS:     if (elt[0] != null) {
    // JTS:       buf.append("A:");
    // JTS:       buf.append(elt[0].toString());
    // JTS:     }
    // JTS:     if (elt[1] != null) {
    // JTS:       buf.append(" B:");
    // JTS:       buf.append(elt[1].toString());
    // JTS:     }
    // JTS:     return buf.toString();
    // JTS:   }
    // JTS: }
}
