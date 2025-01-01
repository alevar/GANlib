extern crate bio;
extern crate petgraph;

// pub mod graph;
pub mod treader;
pub mod factory;
pub mod object;
pub mod transcript;
pub mod exon;
pub mod bundle;
pub mod utils;

pub mod prelude {
    // pub use crate::graph::{GffObject};
    pub use crate::treader::TReader;
    pub use crate::factory::GffObjectFactory;
    pub use crate::object::GffObjectT;
    pub use crate::transcript::Transcript;
    pub use crate::exon::Exon;
    pub use crate::bundle::Bundle;
}

pub use prelude::*;

#[cfg(test)]
mod tests {
}