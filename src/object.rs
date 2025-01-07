// defines traits and srtuct for the general purpose GffObject

use std::convert::TryFrom;
use std::error::Error;

use bio::utils::Interval;
use bio::data_structures::interval_tree::EntryT;

use std::collections::HashMap;
use std::cmp::Ordering;

use crate::utils::*;

pub trait GffObjectT: EntryT<N = usize> + std::fmt::Debug {
    fn seqid(&self) -> &str;
    fn strand(&self) -> char;
    fn get_type(&self) -> Types;
    fn source(&self) -> &str;

    fn score(&self) -> Option<f32> {
        None
    }
    fn phase(&self) -> Option<u32> {
        None
    }

    fn get_attr(&self, key: &str) -> Option<&String> {
        self.get_attrs().get(key)
    }

    fn set_attr(&mut self, key: &str, value: String);

    fn get_attrs(&self) -> &HashMap<String, String>;

    fn bed(&self) -> String;
    fn gtf(&self) -> String;
    fn gff(&self) -> String;

    // len and overlaps/contains are now automatically available from EntryT's interval
    fn len(&self) -> usize {
        (self.interval().end - self.interval().start) + 1
    }

    fn children(&self) -> &[usize];

    fn set_type(&mut self, gtype: Types);

    fn overlaps(&self, other: &dyn GffObjectT) -> bool {
        self.seqid() == other.seqid()
            && self.interval().start <= other.interval().end
            && self.interval().end >= other.interval().start
    }

    fn contains(&self, other: &dyn GffObjectT) -> bool {
        self.seqid() == other.seqid()
            && self.interval().start <= other.interval().start
            && self.interval().end >= other.interval().end
    }
}

// implement a generic object type which can then be specialized into anything
#[derive(Clone, Debug)]
pub struct GffObject {
    pub seqid: String,
    pub strand: char,
    pub source: String,
    pub g_type: Types,
    pub attrs: HashMap<String, String>,
    extra_attrs: HashMap<String,String>, // extra attributes that are not part of the GFF/GTF 9th column
    
    pub id: Option<usize>,
    pub interval: Interval<usize>,
    pub children: Vec::<usize>,
    pub parent: Option<usize>,
}

impl Default for GffObject {
    fn default() -> Self {
        GffObject {
            seqid: String::new(),
            source: String::from("GANLIB"),
            g_type: Types::Unknown,
            strand: '.',
            attrs: HashMap::new(),
            extra_attrs: HashMap::new(),

            id: None,
            interval: Interval::new(0..0).unwrap(),
            children: Vec::new(),
            parent: None,
        }
    }
}

impl EntryT for GffObject {
    type N = usize;

    fn interval(&self) -> &Interval<Self::N> {
        &self.interval
    }
}

impl<'a> EntryT for &'a GffObject {
    type N = usize;

    fn interval(&self) -> &Interval<Self::N> {
        &self.interval
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
            obj.interval = Interval::new(lcs[3].parse::<usize>().unwrap()..lcs[4].parse::<usize>().unwrap())
                    .unwrap();
            obj.strand = lcs[6].chars().next().unwrap();

            obj.g_type = match lcs[2].to_lowercase().as_str() {
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

            obj.attrs = extract_attributes(lcs[8]);

            // add raw source information to the attributes just in case
            obj.extra_attrs = HashMap::new();
            obj.extra_attrs.insert("record_source".to_string(), lcs[2].to_string());

            return Ok(obj);
        }
    }
}

impl PartialEq<GffObject> for GffObject {
    fn eq(&self, other: &GffObject) -> bool {
        self.strand == other.strand()
            && self.seqid == other.seqid()
            && self.interval == *other.interval()
            && self.source == other.source()
            && self.g_type == other.get_type()
            && self.attrs == *other.get_attrs()
    }
}

impl Eq for GffObject {}

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

impl GffObjectT for GffObject {
    fn seqid(&self) -> &str {
        &self.seqid
    }
    fn strand(&self) -> char {
        self.strand
    }
    fn get_type(&self) -> Types {
        self.g_type.clone()
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
                self.g_type,
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
                self.g_type,
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
    fn set_attr(&mut self, key: &str, value: String) {
        self.attrs.insert(key.to_string(), value);
    }

    fn children(&self) -> &[usize] {
        &self.children
    }

    fn set_type(&mut self, gtype: Types) {
        self.g_type = gtype;
    }
}

impl GffObject {
    pub fn new(line: &str) -> Result<GffObject, Box<dyn Error>> {
        GffObject::try_from(line)
    }
}

#[cfg(test)]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_object_creation() {
        let line = "chr1\ttest\tgene\t100\t200\t.\t+\t.\tgene_id \"test\"; gene_name \"test\";";
        let obj = GffObject::new(line).unwrap();
        assert_eq!(obj.seqid, "chr1");
        assert_eq!(obj.source, "test");
        assert_eq!(obj.g_type, Types::Gene);
        assert_eq!(obj.interval.start, 100);
        assert_eq!(obj.interval.end, 200);
        assert_eq!(obj.strand, '+');
        assert_eq!(obj.attrs.len(), 2);
        assert_eq!(obj.attrs.get("gene_id").unwrap(), "test");
        assert_eq!(obj.attrs.get("gene_name").unwrap(), "test");
    }
}