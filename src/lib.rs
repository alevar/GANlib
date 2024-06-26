#![forbid(unsafe_code)]

extern crate bio;

mod object;
pub mod transcript;
pub mod exon;

pub mod treader;

pub mod bundle;
pub mod gene;
pub mod transcriptome;

pub(crate) mod factory;

pub mod prelude {
    pub use crate::transcript::Transcript;
    pub use crate::exon::Exon;

    pub use crate::treader::TReader;

    pub use crate::object::{GffObject,Types};
    use crate::object::GffObjectT;

    pub use crate::gene::Gene;
    pub use crate::transcriptome::Transcriptome;
}

pub use prelude::*;
pub use bundle::Bundle;
pub use factory::GffObjectFactory;

#[cfg(test)]
mod tests {
    use super::*;

    // Unit tests
    #[test]
    fn test() {
    }
}