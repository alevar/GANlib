pub mod utils;
pub mod object;
pub mod group;
pub mod transcript;
pub mod treader;

pub mod prelude {
    pub use crate::object::GffObjectT;
    pub use crate::group::{GffObjectGroupT, Transcriptome};
    pub use crate::transcript::TranscriptRef;
    pub use crate::utils::*;
    pub use crate::treader::TReader;
}

pub use prelude::*;

#[cfg(test)]
mod tests {
}