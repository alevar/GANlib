extern crate bio;
extern crate petgraph;

// pub mod treader;
// pub mod graph;
pub mod object;
pub mod utils;

pub mod prelude {
    // pub use crate::treader::TReader;
    // pub use crate::graph::{GffObject};
    pub use crate::object::{GffObjectT};
}

pub use prelude::*;

#[cfg(test)]
mod tests {
    use super::*;

    // Unit tests
    #[test]
    fn test() {
    }
}