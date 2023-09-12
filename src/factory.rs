use crate::object::{GffObject, GffObjectT, Types};
use crate::transcript::Transcript;
use crate::exon::Exon;
use crate::cds::CDS;

pub struct GTFObjectFactory;

impl GTFObjectFactory {
    fn create(line: &str) -> impl GffObjectT {
        let mut obj = GffObject::new(line).unwrap();
        match obj.get_type() {
            Types::Transcript => Transcript::default(),//obj.to_transcript(),
            Types::Exon => Transcript::default(),//.to_exon(),
            Types::CDS => Transcript::default(),//CDS::new(obj),
            _ => Transcript::default(),
        }
    }
}