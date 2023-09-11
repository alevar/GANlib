
use std::collections::HashMap;
use bio::data_structures::interval_tree::IntervalTree;
use bio::utils::Interval;
use std::cmp::Ordering;
use std::fmt::{Formatter,Display};


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
impl Display for Types{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Types::Gene => write!(f, "gene"),
            Types::Transcript => write!(f, "transcript"),
            Types::Exon => write!(f, "exon"),
            Types::CDS => write!(f, "CDS"),
            Types::UTR => write!(f, "UTR"),
            Types::Intron => write!(f, "intron"),
            Types::Intergenic => write!(f, "intergenic"),
            Types::Unknown => write!(f, "unknown"),
        }
    }
}

pub(crate) trait ObjectT: Clone + Ord + PartialOrd + Eq + PartialEq + Default {
    fn new(line: &str) -> Option<Self>; // return None if line is invalid
    fn add_line(&mut self, line: &str) -> Option<bool>; // returns None if line is invalid, Some(true) if line is valid and object was empty and info was stored
    fn start(&self) -> u32{
        self.interval().start
    }
    fn end(&self) -> u32{
        self.interval().end
    }
    fn get_type(&self) -> Types;
    fn source(&self) -> &str;
    fn score(&self) -> Option<f32>{
        None
    }
    fn phase(&self) -> Option<u32>{
        None
    }
    // fn id(&self) -> Option<&str>;
    // fn attrs(&self) -> &HashMap<String, String>;
    // fn set_attr(&mut self, key: &str, value: String);
    fn bed(&self) -> String;
    fn gtf(&self) -> String;
    fn gff(&self) -> String;
    // fn equals<T: ObjectT>(&self, other: &T) -> bool;
    // fn to_transcript(&self) -> Transcript;
    // fn to_exon(&self) -> Exon;
    fn interval(&self) -> &Interval<u32>;
}


// implement a generic object type which can then be specialized into anything
#[derive(Debug, Clone)]
pub(crate) struct Object {
    seqid: String,
    strand: char,
    interval: Interval<u32>,
    source: String,
    obj_type: Types,
    attrs: HashMap<String, String>,
}

impl ObjectT for Object {
    fn new(line: &str) -> Option<Self> {
        let mut obj = Object::default();
        obj.add_line(line)?;
        Some(obj)
    }
    fn add_line(&mut self, line: &str) -> Option<bool> {
        let lcs: Vec<&str> = line.split('\t').collect();
        if lcs.len() != 9 {
            return None;
        }
        self.interval = Interval::new(lcs[3].parse::<u32>().ok()?..lcs[4].parse::<u32>().ok()?).ok()?;
        self.source = lcs[1].to_string();
        self.obj_type = match lcs[2] {
            "gene" => Types::Gene,
            "transcript" => Types::Transcript,
            "exon" => Types::Exon,
            "CDS" => Types::CDS,
            "UTR" => Types::UTR,
            "intron" => Types::Intron,
            "intergenic" => Types::Intergenic,
            _ => Types::Unknown,
        };
        let mut attrs = HashMap::new();
        for attr in lcs[8].split(';') {
            let kv: Vec<&str> = attr.split('=').collect();
            if kv.len() != 2 {
                continue;
            }
            attrs.insert(kv[0].to_string(), kv[1].to_string());
        }
        Some(true)
    }
    fn interval(&self) -> &Interval<u32> {
        &self.interval
    }
    fn get_type(&self) -> Types {
        self.obj_type.clone()
    }
    fn source(&self) -> &str {
        &self.source
    }
    fn bed(&self) -> String {
        format!("{}\t{}\t{}\t{}\t{}\t{}",
                self.seqid,
                self.interval.start,
                self.interval.end,
                self.source,
                self.score().unwrap_or(0.0),
                self.strand)
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

impl Default for Object {
    fn default() -> Self {
        Object {
            seqid: String::new(),
            interval: Interval::new(0..0).unwrap(),
            source: String::from("GANLIB"),
            obj_type: Types::Unknown,
            attrs: HashMap::new(),
            strand: '.',
        }
    }
}

impl Eq for Object {}

impl PartialEq for Object {
    fn eq(&self, other: &Self) -> bool {
        self.strand == other.strand
            && self.seqid == other.seqid
            && self.interval == other.interval
            && self.source == other.source
            && self.obj_type == other.obj_type
            && self.attrs == other.attrs
    }
}

impl Ord for Object {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.interval.start.cmp(&other.interval.start) { // TODO: need to add seqid and strand
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

