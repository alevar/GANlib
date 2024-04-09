use std::collections::HashMap;
use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use bio::utils::Interval;
use std::cmp::Ordering;
use std::fmt::{Formatter,Display};
use std::error::Error;


use crate::bundle::GffObjectGroupT;
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

pub trait GffObjectT {
    // fn new(line: &str) -> Result<Self,Box<dyn Error>>; // return None if line is invalid
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
    fn len(&self) -> u32{
        (self.interval().end - self.interval().start)+1
    }
    // fn overlaps(&self, other: &Self) -> bool{
    //     self.seqid() == other.seqid() && self.interval().overlaps(other.interval())
    // }
    // fn contains(&self, other: &Self) -> bool{
    //     self.seqid() == other.seqid() && self.interval().contains(other.interval())
    // }
    
}


// implement a generic object type which can then be specialized into anything
#[derive(Debug, Clone)]
pub struct GffObject {
    pub seqid: String,
    pub strand: char,
    pub interval: Interval<u32>,
    pub source: String,
    pub obj_type: Types,
    pub attrs: HashMap<String, String>,
    extra_attrs: HashMap<String,String>, // extra attributes that are not part of the GFF/GTF 9th column
    children: ArrayBackedIntervalTree<u32,GffObject>, // TODO: add children - any object should be able to inherit children when converted from something like transcript or gene, or anyhting else. It might not be able to do anything with the children, but should be able to store them
}

impl GffObject{
    pub fn new(line: &str) -> Result<GffObject,Box<dyn Error>> {
        GffObject::try_from(line)
    }
}

impl GffObjectT for GffObject {
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

impl GffObjectGroupT for GffObject {
    fn add(&mut self, obj: Self) {
        self.children.insert(obj.interval().clone(), obj.clone());
    }
    fn num_elements(&self) -> usize {
        self.children.len()
    }
    fn iter(&self) -> Box<dyn Iterator<Item = &Self>> {
        Box::new(self.children.into_iter())
    }
}

impl From<Exon> for GffObject {
    fn from(exon: Exon) -> GffObject {
        let mut obj = GffObject::default();
        obj.seqid = exon.seqid().to_string();
        obj.strand = exon.strand();
        obj.interval = exon.interval().clone();
        obj.source = exon.source().to_string();
        obj.obj_type = exon.get_type();
        obj.attrs = exon.get_attrs().clone();
        obj
    }
}
impl From<Transcript> for GffObject {
    fn from(tx: Transcript) -> GffObject {
        let mut obj = GffObject::default();
        obj.seqid = tx.seqid().to_string();
        obj.strand = tx.strand();
        obj.interval = tx.interval().clone();
        obj.source = tx.source().to_string();
        obj.obj_type = tx.get_type();
        obj.attrs = tx.get_attrs().clone();
        for exon in tx.exons(){
            obj.children.insert(exon.interval().clone(),GffObject::from(exon.data().clone()));
        }
        obj
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
            children: ArrayBackedIntervalTree::new(),
        }
    }
}

// implementation of from for GffObject conversion from string
impl TryFrom<&str> for GffObject {
    type Error = Box<dyn Error>;
    fn try_from(line: &str) -> Result<Self,Self::Error> {
        // parse line (gtf or gff)
        let mut obj = GffObject::default();
        
        let lcs: Vec<&str> = line.split('\t').collect();
        if lcs.len() != 9 {
            Err(format!("Invalid number of columns in GFF/GTF line: {}", line).into())
        }
        else{
            obj.seqid = lcs[0].to_string();
            obj.source = lcs[1].to_string();
            obj.interval = Interval::new(lcs[3].parse::<u32>().unwrap()..lcs[4].parse::<u32>().unwrap()).unwrap();
            obj.strand = lcs[6].chars().next().unwrap();

            obj.attrs = HashMap::new();
            for attr in lcs[8].split(';') {
                // split on either ' ' or '='
                let mut kv = attr.trim().split(|c| c==' ' || c=='=').collect::<Vec<&str>>();
                if kv.len() != 2 {
                    continue;
                }
                obj.attrs.insert(kv[0].to_lowercase().to_string(), kv[1].trim_matches('"').to_string());
            }

            obj.obj_type = match lcs[2].to_lowercase().as_str() {
                "gene" => Types::Gene,
                "transcript" => Types::Transcript,
                "mrna" => Types::Transcript,
                "exon" => Types::Exon,
                "cds" => Types::CDS,
                "utr" => Types::UTR,
                "intron" => Types::Intron,
                "intergenic" => Types::Intergenic,
                _ => Types::Unknown,
            };

            // process some type-specific policies
            match obj.obj_type {
                // if transcript type - must have transcript_id and gene_id
                Types::Transcript => {
                    // rename potential ID attribute to transcript_id
                    if let Some(id) = obj.attrs.remove("id") {
                        obj.attrs.insert("transcript_id".to_string(), id);
                    }
                    if let Some(id) = obj.attrs.remove("parent") {
                        obj.attrs.insert("gene_id".to_string(), id);
                    }
                },
                Types::Exon => {
                    if let Some(id) = obj.attrs.remove("id") {
                        obj.attrs.insert("exon_id".to_string(), id);
                    }
                    if let Some(id) = obj.attrs.remove("parent") {
                        obj.attrs.insert("transcript_id".to_string(), id);
                    }
                },
                Types::CDS => {
                    if let Some(id) = obj.attrs.remove("id") {
                        obj.attrs.insert("cds_id".to_string(), id);
                    }
                    if let Some(id) = obj.attrs.remove("parent") {
                        obj.attrs.insert("transcript_id".to_string(), id);
                    }
                },
                Types::Gene => {
                    if let Some(id) = obj.attrs.remove("id") {
                        obj.attrs.insert("gene_id".to_string(), id);
                    }
                    if let Some(id) = obj.attrs.remove("parent") {
                        obj.attrs.insert("gene_parent_id".to_string(), id);
                    }
                },
                Types::Intron => {
                    if let Some(id) = obj.attrs.remove("id") {
                        obj.attrs.insert("intron_id".to_string(), id);
                    }
                    if let Some(id) = obj.attrs.remove("parent") {
                        obj.attrs.insert("transcript_id".to_string(), id);
                    }
                },
                Types::UTR => {
                    if let Some(id) = obj.attrs.remove("id") {
                        obj.attrs.insert("utr_id".to_string(), id);
                    }
                    if let Some(id) = obj.attrs.remove("parent") {
                        obj.attrs.insert("transcript_id".to_string(), id);
                    }
                },
                _ => {},
            }
            // add raw source information to the attributes just in case
            obj.extra_attrs = HashMap::new();
            obj.extra_attrs.insert("record_source".to_string(), lcs[2].to_string());

            return Ok(obj);
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
        // compare based on seqid, strand and interval
        match self.seqid.cmp(&other.seqid) {
            Ordering::Equal => match self.strand.cmp(&other.strand) {
                Ordering::Equal => match self.interval.start.cmp(&other.interval.start) {
                    Ordering::Equal => self.interval.end.cmp(&other.interval.end),
                    other => other,
                },
                other => other,
            },
            other => other,
        }
    }
}

impl PartialOrd for GffObject {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // test construction of default empty object and populating it with data after construction
    #[test]
    fn test_gffobject_from_default() {
        let mut obj = GffObject::default();
        assert_eq!(obj.seqid, "");
        assert_eq!(obj.source, "GANLIB");
        assert_eq!(obj.interval, Interval::new(0..0).unwrap());
        assert_eq!(obj.strand, '.');
        assert_eq!(obj.attrs.len(), 0);
        assert_eq!(obj.obj_type, Types::Unknown);
        assert_eq!(obj.extra_attrs.len(), 0);
        obj.seqid = "chr1".to_string();
        obj.source = "test".to_string();
        obj.interval = Interval::new(1..10).unwrap();
        obj.strand = '+';
        obj.attrs.insert("test".to_string(), "test".to_string());
        obj.obj_type = Types::Gene;
        obj.extra_attrs.insert("test".to_string(), "test".to_string());
        assert_eq!(obj.seqid, "chr1");
        assert_eq!(obj.source, "test");
        assert_eq!(obj.interval, Interval::new(1..10).unwrap());
        assert_eq!(obj.strand, '+');
        assert_eq!(obj.attrs.len(), 1);
        assert_eq!(obj.attrs.get("test").unwrap(), "test");
        assert_eq!(obj.obj_type, Types::Gene);
        assert_eq!(obj.extra_attrs.len(), 1);
        assert_eq!(obj.extra_attrs.get("test").unwrap(), "test");
    }

    #[test]
    fn test_gffobject() {
        let line = "chr1\t.\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\"; gene_type \"transcribed_unprocessed_pseudogene\"; gene_name \"DDX11L1\"; level 2; havana_gene \"OTTHUMG00000000961.2\"; remap_status \"full_contig\";";
        let obj = GffObject::new(line).unwrap();
        assert_eq!(obj.seqid, "chr1");
        assert_eq!(obj.source, ".");
        assert_eq!(obj.interval, Interval::new(11869..14409).unwrap());
        assert_eq!(obj.strand, '+');
        assert_eq!(obj.attrs.len(), 6);
        assert_eq!(obj.attrs.get("gene_id").unwrap(), "ENSG00000223972.5");
        assert_eq!(obj.attrs.get("gene_type").unwrap(), "transcribed_unprocessed_pseudogene");
        assert_eq!(obj.attrs.get("gene_name").unwrap(), "DDX11L1");
        assert_eq!(obj.attrs.get("level").unwrap(), "2");
        assert_eq!(obj.attrs.get("havana_gene").unwrap(), "OTTHUMG00000000961.2");
        assert_eq!(obj.attrs.get("remap_status").unwrap(), "full_contig");
        assert_eq!(obj.obj_type, Types::Gene);
        assert_eq!(obj.extra_attrs.get("record_source").unwrap(), "gene");
    }
    
    // test equality of two GffObjects
    #[test]
    fn test_gffobject_eq() {
        // test equality
        let line = "chr1\t.\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\"; gene_type \"transcribed_unprocessed_pseudogene\"; gene_name \"DDX11L1\"; level 2; havana_gene \"OTTHUMG00000000961.2\"; remap_status \"full_contig\";";
        let obj1 = GffObject::new(line).unwrap();
        let obj2 = GffObject::new(line).unwrap();
        assert_eq!(obj1, obj2);
        // test inequality
        let line = "chr1\t.\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\"; gene_type \"transcribed_unprocessed_pseudogene\"; gene_name \"DDX11L1\"; level 2; havana_gene \"OTTHUMG00000000961.2\";";
        let obj3 = GffObject::new(line).unwrap();
        assert_ne!(obj1, obj3);
    }
    
    // test ordering of two GffObjects
    #[test]
    fn test_gffobject_ord() {
        // test equality, and ordering based on the seqid
        let line = "chr1\t.\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\"; gene_type \"transcribed_unprocessed_pseudogene\"; gene_name \"DDX11L1\"; level 2; havana_gene \"OTTHUMG00000000961.2\"; remap_status \"full_contig\";";
        let obj1 = GffObject::new(line).unwrap();
        let line = "chr2\t.\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\"; gene_type \"transcribed_unprocessed_pseudogene\"; gene_name \"DDX11L1\"; level 2; havana_gene \"OTTHUMG00000000961.2\";";
        let obj2 = GffObject::new(line).unwrap();
        assert!(obj1 < obj2);

        // test based on strand
        let line = "chr1\t.\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\"; gene_type \"transcribed_unprocessed_pseudogene\"; gene_name \"DDX11L1\"; level 2; havana_gene \"OTTHUMG00000000961.2\"; remap_status \"full_contig\";";
        let obj1 = GffObject::new(line).unwrap();
        let line = "chr1\t.\tgene\t11869\t14409\t.\t-\t.\tgene_id \"ENSG00000223972.5\"; gene_type \"transcribed_unprocessed_pseudogene\"; gene_name \"DDX11L1\"; level 2; havana_gene \"OTTHUMG00000000961.2\";";
        let obj2 = GffObject::new(line).unwrap();
        assert!(obj1 < obj2);

        // test based on interval
        let line = "chr1\t.\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\"; gene_type \"transcribed_unprocessed_pseudogene\"; gene_name \"DDX11L1\"; level 2; havana_gene \"OTTHUMG00000000961.2\"; remap_status \"full_contig\";";
        let obj1 = GffObject::new(line).unwrap();
        let line = "chr1\t.\tgene\t11869\t14410\t.\t+\t.\tgene_id \"ENSG00000223972.5\"; gene_type \"transcribed_unprocessed_pseudogene\"; gene_name \"DDX11L1\";";
        let obj2 = GffObject::new(line).unwrap();
        assert!(obj1 < obj2);

        // test seqid first and then strand
        let line = "chr1\t.\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\";";
        let obj1 = GffObject::new(line).unwrap();
        let line = "chr1\t.\tgene\t11869\t14409\t.\t-\t.\tgene_id \"ENSG00000223972.5\";";
        let obj2 = GffObject::new(line).unwrap();
        assert!(obj1 < obj2);
        
        // test seqid first, strand second and interval last
        let line = "chr1\t.\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\";";
        let obj1 = GffObject::new(line).unwrap();
        let line = "chr1\t.\tgene\t11869\t14410\t.\t+\t.\tgene_id \"ENSG00000223972.5\";";
        let obj2 = GffObject::new(line).unwrap();
        assert!(obj1 < obj2);
    }

    // test conversion to transcript
    #[test]
    fn test_gffobject_to_transcript() {
        let line = "chr1\t.\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\";";
        let obj = GffObject::new(line).unwrap();
        let tx = Transcript::from(obj);
        assert_eq!(tx.seqid(), "chr1");
        assert_eq!(tx.strand(), '+');
        assert_eq!(tx.interval(), &Interval::new(11869..14409).unwrap());
        assert_eq!(tx.get_type(), Types::Transcript);
        assert_eq!(tx.get_attr("gene_id").unwrap(), "ENSG00000223972.5");
        assert_eq!(tx.get_attr("transcript_id"), None);
        assert_eq!(tx.get_attr("gene_type"), None);
        assert_eq!(tx.get_attr("gene_name"), None);
        assert_eq!(tx.get_attr("level"), None);
        assert_eq!(tx.get_attr("havana_gene"), None);
        assert_eq!(tx.get_attr("remap_status"), None);
    }
}