use crate::object::{GffObject, GffObjectT, Types};
use crate::transcript::Transcript;
use crate::exon::Exon;
use bio::utils::Interval;
use std::error::Error;

pub struct GffObjectFactory;

impl Default for GffObjectFactory {
    fn default() -> Self {
        GffObjectFactory
    }
}

impl GffObjectFactory {
    pub(crate) fn create(&self, line: &str) -> Result<Box<dyn GffObjectT>,Box<dyn Error>> {
        let obj = match GffObject::new(line) {
            Ok(obj) => obj,
            Err(e) => return Err(e),
        };
        match obj.get_type() {
            Types::Transcript => Ok(Box::new(Transcript::from(obj))),
            Types::Exon => Ok(Box::new(Exon::from(obj))),
            Types::CDS => {
                            let mut exon = Exon::from(obj);
                            exon.set_cds(Interval::new(exon.interval().start..exon.interval().end).unwrap());
                            Ok(Box::new(exon))                        
                          },
            _ => Ok(Box::new(obj)),
        }
    }
}