use std::{collections::HashMap, cmp::Ordering};
use bio::io::gff;
use bio::utils::Interval;
use std::error::Error;
use std::ptr::NonNull;

use crate::object::*;
use crate::utils::*;
use crate::exon::Exon;

use bio::data_structures::interval_tree::{ArrayBackedIntervalTree, EntryT};

#[derive(Clone, Debug)]
pub struct Transcript {
    base: GffObject,
    transcript_id: Option<String>,
    gene_id: Option<String>,
    exons: ArrayBackedIntervalTree<Exon>,
}

impl Transcript {
    // Constructor to create an Transcript from a GffObject
    pub fn from(gff_obj: GffObject) -> Self {
        let mut new_transcript = Self {
            base: gff_obj.clone(),
            transcript_id: None,
            gene_id: None,
            exons: ArrayBackedIntervalTree::new(),
        };

        new_transcript
    }
}

impl GffObjectT for Transcript {
    fn new(line: &str) -> Result<Transcript, Box<dyn Error>> {
        Transcript::try_from(line)
    }
    fn parent(&self) -> Option<&dyn GffObjectT> {
        self.base.parent.and_then(|parent_ptr| unsafe { Some(parent_ptr.as_ref() as &dyn GffObjectT) })
    }
    fn seqid(&self) -> &str {
        &self.base.seqid
    }
    fn strand(&self) -> char {
        self.base.strand
    }
    fn get_type(&self) -> Types {
        Types::Transcript
    }
    fn source(&self) -> &str {
        &self.base.source
    }
    fn bed(&self) -> String {
        format!("{}\t{}\t{}\t{}\t{}\t{}",
                self.base.seqid,
                self.base.interval.start,
                self.base.interval.end,
                self.base.source,
                self.base.score().unwrap_or(0.0),
                self.base.strand)
    }
    fn gtf(&self) -> String {
        format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.base.seqid,
                self.base.source,
                Types::Transcript,
                self.base.interval.start,
                self.base.interval.end,
                self.score().unwrap_or(0.0),
                self.base.strand,
                self.phase().unwrap_or(0),
                self.base.attrs.iter().map(|(k,v)| format!("{} \"{}\";", k, v)).collect::<Vec<String>>().join(" "))
    }
    fn gff(&self) -> String {
        format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.base.seqid,
                self.base.source,
                Types::Transcript,
                self.base.interval.start,
                self.base.interval.end,
                self.score().unwrap_or(0.0),
                self.base.strand,
                self.phase().unwrap_or(0),
                self.base.attrs.iter().map(|(k,v)| format!("{}={};", k, v)).collect::<Vec<String>>().join(" "))
    }
    fn get_attrs(&self) -> &HashMap<String, String> {
        &self.base.attrs
    }
    fn get_attrs_mut(&mut self) -> &mut HashMap<String, String> {
        &mut self.base.attrs
    }
}

impl TryFrom<&str> for Transcript {
    type Error = Box<dyn Error>;
    fn try_from(line: &str) -> Result<Self,Self::Error> {
        let mut obj = GffObject::try_from(line)?;
        return Ok(Transcript::from(obj));
    }
}

impl EntryT for Transcript {
    type N = usize;

    fn interval(&self) -> &Interval<Self::N> {
        &self.base.interval
    }
}

impl<'a> EntryT for &'a Transcript {
    type N = usize;

    fn interval(&self) -> &Interval<Self::N> {
        &self.base.interval
    }
}