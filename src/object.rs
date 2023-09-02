
use std::collections::HashMap;
use bio::data_structures::interval_tree::IntervalTree;
use bio::utils::Interval;
use std::cmp::Ordering;


use crate::transcript::Transcript;
use crate::exon::Exon;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Types {
    Gene,
    Transcript,
    Exon,
    CDS,
    UTR,
    Intron,
    Intergenic,
    Unknown,
}

trait ObjectT: Clone + Ord {
    fn get_start(&self) -> u32{
        self.get_interval().start
    }
    fn get_end(&self) -> u32{
        self.get_interval().end
    }
    fn get_type(&self) -> Types;
    fn get_source(&self) -> &str;
    fn get_score(&self) -> Option<f32>{
        None
    }
    fn get_phase(&self) -> Option<u32>{
        None
    }
    fn get_id(&self) -> Option<&str>;
    fn get_attrs(&self) -> &HashMap<String, String>;
    fn set_attr(&mut self, key: &str, value: String);
    fn to_bed_record(&self) -> String;
    fn to_gtf(&self) -> String;
    fn to_gff(&self) -> String;
    fn equals<T: ObjectT>(&self, other: &T) -> bool;
    fn to_transcript(&self) -> Transcript;
    fn to_exon_opt(&self) -> Option<Exon>;
    fn get_interval(&self) -> &Interval<u32>;
}


// implement a generic object type which can then be specialized into anything
#[derive(Debug, Clone)]
pub(crate) struct Object {
    interval: Interval<u32>,
    source: String,
    obj_type: String,
    attrs: HashMap<String, String>,
}

impl Eq for Object {}

impl PartialEq for Object {
    fn eq(&self, other: &Self) -> bool {
        self.interval == other.interval
            && self.source == other.source
            && self.obj_type == other.obj_type
            && self.attrs == other.attrs
    }
}

impl Ord for Object {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.interval.start.cmp(&other.interval.start) {
            Ordering::Equal => self.interval.end.cmp(&other.interval.end),
            other => other,
        }
    }
}

impl PartialOrd for Object {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}