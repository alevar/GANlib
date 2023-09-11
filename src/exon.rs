use std::{collections::HashMap, cmp::Ordering};
use bio::utils::Interval;
use std::ops::Deref;

use crate::object::*
;

#[derive(Clone, Debug)]
pub struct Exon {
    seqid: String,
    strand: char,
    interval: Interval<u32>,
    source: String,
    tid: Option<String>,
    gid: Option<String>,
    attrs: HashMap<String, String>,
    obj_type: Types,
}

impl Default for Exon {
    fn default() -> Self {
        Exon {
            seqid: String::new(),
            strand: '.',
            source: String::from("GANLIB"),
            interval: Interval::new(0..0).unwrap(),
            tid: None,
            gid: None,
            attrs: HashMap::new(),
            obj_type: Types::Exon,
        }
    }
}

impl Exon {
    pub fn new() -> Self {
        Exon {
            seqid: String::new(),
            strand: '.',
            interval: Interval::default(),
            source: String::from("GANLIB"),
            tid: None,
            gid: None,
            attrs: HashMap::new(),
            obj_type: Types::Exon,
        }
    }

    pub fn set_tid(&mut self, tid: String) {
        self.tid = Some(tid);
    }

    pub fn set_gid(&mut self, gid: String) {
        self.gid = Some(gid);
    }

    pub fn get_tid(&self) -> Option<&str> {
        self.tid.as_deref()
    }

    pub fn get_gid(&self) -> Option<&str> {
        self.gid.as_deref()
    }
}

impl ObjectT for Exon {
    fn new(line: &str) -> Option<Self> {
        let mut e = Exon::default();
        e.add_line(line);
        Some(e)
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
        format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.seqid,
                self.source,
                self.obj_type,
                self.interval.start,
                self.interval.end,
                self.score().unwrap_or(0.0),
                self.strand,
                self.phase().unwrap_or(0),
                self.attrs.iter().map(|(k,v)| format!("{} \"{}\";", k, v)).collect::<Vec<String>>().join(" "))
    }
    fn gff(&self) -> String {
        format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                        self.seqid,
                                        self.source,
                                        self.obj_type,
                                        self.interval.start,
                                        self.interval.end,
                                        self.score().unwrap_or(0.0),
                                        self.strand,
                                        self.phase().unwrap_or(0),
                                        self.attrs.iter().map(|(k,v)| format!("{}={};", k, v)).collect::<Vec<String>>().join(" "))
    }
}

impl Eq for Exon {}

impl PartialEq for Exon {
    fn eq(&self, other: &Self) -> bool {
        self.strand == other.strand
            && self.seqid == other.seqid
            && self.interval == other.interval
            && self.source == other.source
            && self.obj_type == other.obj_type
            && self.attrs == other.attrs
    }
}

impl Ord for Exon {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.seqid != other.seqid {
            return self.seqid.cmp(&other.seqid);
        } else if self.strand != other.strand {
            return self.strand.cmp(&other.strand);
        } else {
            match self.interval().end.cmp(&other.interval().end) {
                Ordering::Equal => self.interval().end.cmp(&other.interval().end),
                other => other,
            }
        }

    }
}

impl PartialOrd for Exon {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}