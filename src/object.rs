
use std::collections::HashMap;
use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
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

pub(crate) trait GffObjectT: Clone + Ord + PartialOrd + Eq + PartialEq + Default {
    fn new(line: &str) -> Option<Self>; // return None if line is invalid
    fn start(&self) -> u32{
        self.interval().start
    }
    fn end(&self) -> u32{
        self.interval().end
    }
    fn seqid(&self) -> &str;
    fn strand(&self) -> char;
    fn get_type(&self) -> Types;
    fn source(&self) -> &str;
    fn score(&self) -> Option<f32>{
        None
    }
    fn phase(&self) -> Option<u32>{
        None
    }
    fn get_attr(&self, key: &str) -> Option<&String> {
        self.get_attrs().get(key)
    }
    fn set_attr(&mut self, key: &str, value: String) {
        self.get_attrs_mut().insert(key.to_string(), value);
    }
    fn get_attrs(&self) -> &HashMap<String, String>{
        unimplemented!()
    }
    fn get_attrs_mut(&mut self) -> &mut HashMap<String, String>{
        unimplemented!()
    }
    fn bed(&self) -> String;
    fn gtf(&self) -> String;
    fn gff(&self) -> String;
    fn interval(&self) -> &Interval<u32>;
}


// implement a generic object type which can then be specialized into anything
#[derive(Debug, Clone)]
pub struct GffObject {
    seqid: String,
    strand: char,
    interval: Interval<u32>,
    source: String,
    obj_type: Types,
    attrs: HashMap<String, String>,
    extra_attrs: HashMap<String,String>, // extra attributes that are not part of the GFF/GTF 9th column
}

impl GffObjectT for GffObject {
    fn new(line: &str) -> Option<Self> {
        let mut obj = GffObject::from(line);
        Some(obj)
    }
    fn seqid(&self) -> &str {
        &self.seqid
    }
    fn strand(&self) -> char {
        self.strand
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
    fn get_attrs(&self) -> &HashMap<String, String> {
        &self.attrs
    }
    fn get_attrs_mut(&mut self) -> &mut HashMap<String, String> {
        &mut self.attrs
    }
}

impl From<Exon> for GffObject {
    fn from(exon: Exon) -> GffObject {
        unimplemented!()
    }
}
impl From<Transcript> for GffObject {
    fn from(tx: Transcript) -> GffObject {
        unimplemented!()
    }
}

impl Default for GffObject {
    fn default() -> Self {
        GffObject {
            seqid: String::new(),
            interval: Interval::new(0..0).unwrap(),
            source: String::from("GANLIB"),
            obj_type: Types::Unknown,
            strand: '.',
            attrs: HashMap::new(),
            extra_attrs: HashMap::new(),
        }
    }
}

// implementation of from for GffObject conversion from string
impl From<&str> for GffObject {
    fn from(line: &str) -> Self {
        // parse line (gtf or gff)
        let mut obj = GffObject::default();
        
        let lcs: Vec<&str> = line.split('\t').collect();
        if lcs.len() != 9 {
            panic!("Invalid GFF/GTF line: {}", line);            
        }
        else{
            obj.seqid = lcs[0].to_string();
            obj.source = lcs[1].to_string();
            obj.interval = Interval::new(lcs[3].parse::<u32>().unwrap()..lcs[4].parse::<u32>().unwrap()).unwrap();
            obj.strand = lcs[6].chars().next().unwrap();

            obj.attrs = HashMap::new();
            for attr in lcs[8].split(';') {
                // split on either ' ' or '='
                let mut kv = attr.split(|c| c==' ' || c=='=').collect::<Vec<&str>>();
                if kv.len() != 2 {
                    continue;
                }
                obj.attrs.insert(kv[0].to_string(), kv[1].to_string());
            }

            obj.obj_type = match lcs[2] {
                "gene" => Types::Gene,
                "transcript" => Types::Transcript,
                "exon" => Types::Exon,
                "CDS" => Types::CDS,
                "UTR" => Types::UTR,
                "intron" => Types::Intron,
                "intergenic" => Types::Intergenic,
                _ => Types::Unknown,
            };
            // add raw source information to the attributes just in case
            obj.extra_attrs = HashMap::new();
            obj.extra_attrs.insert("record_source".to_string(), lcs[2].to_string());

            return obj;
        }
    }
}


impl Eq for GffObject {}

impl PartialEq for GffObject {
    fn eq(&self, other: &Self) -> bool {
        self.strand == other.strand
            && self.seqid == other.seqid
            && self.interval == other.interval
            && self.source == other.source
            && self.obj_type == other.obj_type
            && self.attrs == other.attrs
    }
}

impl Ord for GffObject {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.interval.start.cmp(&other.interval.start) { // TODO: need to add seqid and strand
            Ordering::Equal => self.interval.end.cmp(&other.interval.end),
            other => other,
        }
    }
}

impl PartialOrd for GffObject {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

