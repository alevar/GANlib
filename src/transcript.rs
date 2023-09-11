use std::collections::HashMap;
use std::fmt::format;
use bio::utils::Interval;
use bio::data_structures::interval_tree::{IntervalTree, IntervalTreeIterator, ArrayBackedIntervalTree};
use std::cmp::Ordering;
use std::ops::Deref;

use crate::object::{Object, ObjectT, Types};
use crate::exon::Exon;


#[derive(Clone, Debug)]
pub struct Transcript {
    seqid: String,
    strand: char,
    source: String,
    interval: Interval<u32>,
    exons: ArrayBackedIntervalTree<u32,Exon>,
    tid: Option<String>,
    gid: Option<String>,
    attrs: HashMap<String, String>,
    obj_type: Types,
}

impl Default for Transcript {
    fn default() -> Self {
        Transcript {
            seqid: String::new(),
            strand: '.',
            source: String::from("GANLIB"),
            interval: Interval::new(0..0).unwrap(),
            exons: ArrayBackedIntervalTree::new(),
            tid: None,
            gid: None,
            attrs: HashMap::new(),
            obj_type: Types::Transcript,
        }
    }
}

impl ObjectT for Transcript {
    fn new(line: &str) -> Option<Self> {
        let mut t = Transcript::default();
        t.add_line(line);
        Some(t)
    }

    fn add_line(&mut self, line: &str) -> Option<bool> {
        Some(true)
    }

    fn get_type(&self) -> Types {
        Types::Transcript
    }

    fn interval(&self) -> &Interval<u32> {
        &self.interval
    }
    fn source(&self) -> &str {
        &self.source
    }
    fn bed(&self) -> String {
        format!("{}\t{}\t{}\t{}\t{}\t{}",
                self.seqid,
                self.interval.start,
                self.interval.end,
                self.tid.as_ref().unwrap_or(&String::from(".")),
                self.score().unwrap_or(0.0),
                self.strand,
        )
    }
    fn gtf(&self) -> String {
        let mut gtf_str = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.seqid,
                self.source,
                self.obj_type,
                self.interval.start,
                self.interval.end,
                self.score().unwrap_or(0.0),
                self.strand,
                self.phase().unwrap_or(0),
                self.attrs.iter().map(|(k,v)| format!("{} \"{}\";", k, v)).collect::<Vec<String>>().join(" "));
        for exon in &self.exons{
            gtf_str.push_str(&exon.data().gff());
        }
        gtf_str
    }

    fn gff(&self) -> String {
        let mut gff_str = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                        self.seqid,
                                        self.source,
                                        self.obj_type,
                                        self.interval.start,
                                        self.interval.end,
                                        self.score().unwrap_or(0.0),
                                        self.strand,
                                        self.phase().unwrap_or(0),
                                        self.attrs.iter().map(|(k,v)| format!("{}={};", k, v)).collect::<Vec<String>>().join(" "));
        for exon in &self.exons{
            gff_str.push_str(&exon.data().gff());
        }
        gff_str
    }
}

impl Eq for Transcript {}

impl PartialEq for Transcript {
    fn eq(&self, other: &Self) -> bool {
        self.strand == other.strand
            && self.seqid == other.seqid
            && self.exons == other.exons
            && self.source == other.source
            && self.obj_type == other.obj_type
            && self.attrs == other.attrs
    }
}

impl Ord for Transcript {
    fn cmp(&self, other: &Self) -> Ordering {
        self.exons.into_iter().cmp(other.exons.into_iter())
    }
}

impl PartialOrd for Transcript {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}