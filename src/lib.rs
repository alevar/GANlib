#![forbid(unsafe_code)]

extern crate bio;

mod object;
pub mod transcript;
pub mod exon;
pub mod cds;

pub mod treader;

mod txgroup;
pub mod bundle;
pub mod gene;
pub mod transcriptome;

mod prelude {
    pub use crate::transcript::Transcript;
    pub use crate::exon::Exon;
    pub use crate::cds::CDS;

    pub use crate::treader::TReader;

    use crate::object::Object;

    pub use crate::gene::Gene;
    pub use crate::transcriptome::Transcriptome;
}

pub use prelude::*;
pub use txgroup::TXGroup;