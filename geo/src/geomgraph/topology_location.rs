use super::{Location, Position};

use std::fmt;

#[derive(Clone)]
pub(crate) enum TopologyLocation {
    Area {
        on: Option<Location>,
        left: Option<Location>,
        right: Option<Location>,
    },
    Line {
        on: Option<Location>,
    },
}

impl fmt::Debug for TopologyLocation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fn location_to_str(location: &Option<Location>, f: &mut fmt::Formatter) -> fmt::Result {
            match location {
                Some(Location::Interior) => write!(f, "i"),
                Some(Location::Boundary) => write!(f, "b"),
                Some(Location::Exterior) => write!(f, "e"),
                None => write!(f, "_"),
            }
        }
        match self {
            Self::Line { on } => location_to_str(on, f)?,
            Self::Area { on, left, right } => {
                location_to_str(left, f)?;
                location_to_str(on, f)?;
                location_to_str(right, f)?;
            }
        }
        Ok(())
    }
}

impl TopologyLocation {
    pub fn area(on: Location, left: Location, right: Location) -> Self {
        Self::Area {
            on: Some(on),
            left: Some(left),
            right: Some(right),
        }
    }

    pub fn empty_area() -> Self {
        Self::Area {
            on: None,
            left: None,
            right: None,
        }
    }

    pub fn line(on: Location) -> Self {
        Self::Line { on: Some(on) }
    }

    pub fn empty_line() -> Self {
        Self::Line { on: None }
    }

    pub fn get(&self, pos: Position) -> Option<Location> {
        match (pos, self) {
            (Position::On, Self::Line { on }) | (Position::On, Self::Area { on, .. }) => *on,
            (Position::Left, Self::Area { left, .. }) => *left,
            (Position::Right, Self::Area { right, .. }) => *right,
            _ => {
                debug_assert!(false, "does this happen?");
                return None;
            }
        }
    }

    pub fn is_empty(&self) -> bool {
        match self {
            Self::Line { on: None } => true,
            Self::Area {
                on: None,
                left: None,
                right: None,
            } => true,
            _ => false,
        }
    }

    pub fn is_any_empty(&self) -> bool {
        match self {
            Self::Line { on: Some(_) } => false,
            Self::Area {
                on: Some(_),
                left: Some(_),
                right: Some(_),
            } => false,
            _ => true,
        }
    }

    pub fn is_area(&self) -> bool {
        matches!(self, Self::Area { .. })
    }

    pub fn is_line(&self) -> bool {
        matches!(self, Self::Line { .. })
    }

    pub fn flip(&mut self) {
        match self {
            Self::Line { .. } => {}
            Self::Area { left, right, .. } => {
                std::mem::swap(left, right);
            }
        }
    }

    pub fn set_all_locations(&mut self, location: Location) {
        match self {
            Self::Line { on } => {
                *on = Some(location);
            }
            Self::Area { on, left, right } => {
                *on = Some(location);
                *left = Some(location);
                *right = Some(location);
            }
        }
    }

    pub fn set_all_locations_if_empty(&mut self, location: Location) {
        match self {
            Self::Line { on } => {
                if on.is_none() {
                    *on = Some(location);
                }
            }
            Self::Area { on, left, right } => {
                if on.is_none() {
                    *on = Some(location);
                }
                if left.is_none() {
                    *left = Some(location);
                }
                if right.is_none() {
                    *right = Some(location);
                }
            }
        }
    }

    pub fn set_location(&mut self, position: Position, location: Location) {
        match (position, self) {
            (Position::On, Self::Line { on }) => *on = Some(location),
            (_, Self::Line { .. }) => {
                panic!("invalid assignment dimensions for Self::Line")
            }
            (Position::On, Self::Area { on, .. }) => *on = Some(location),
            (Position::Left, Self::Area { left, .. }) => *left = Some(location),
            (Position::Right, Self::Area { right, .. }) => *right = Some(location),
        }
    }

    pub fn set_on_location(&mut self, location: Location) {
        match self {
            Self::Line { on } | Self::Area { on, .. } => {
                *on = Some(location);
            }
        }
    }

    pub fn set_locations(&mut self, new_on: Location, new_left: Location, new_right: Location) {
        match self {
            Self::Line { .. } => {
                error!("invalid assignment dimensions for {:?}", self);
                debug_assert!(false, "invalid assignment dimensions for {:?}", self);
            }
            Self::Area { on, left, right } => {
                *on = Some(new_on);
                *left = Some(new_left);
                *right = Some(new_right);
            }
        }
    }
}
