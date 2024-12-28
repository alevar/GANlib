use std::convert::TryFrom;
use std::error::Error;
use std::ptr::NonNull;

use bio::utils::Interval;
use bio::data_structures::interval_tree::{ArrayBackedIntervalTree, EntryT};

use std::collections::HashMap;
use std::cmp::Ordering;

use crate::utils::*;

pub trait GffObjectT: EntryT<N = usize> {
    // Create a new instance from a GFF line
    fn new(line: &str) -> Result<Self, Box<dyn Error>>
    where
        Self: Sized;

    // Returns the parent as an optional pointer
    fn parent(&self) -> Option<&dyn GffObjectT>;

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

    fn set_attr(&mut self, key: &str, value: String) {
        self.get_attrs_mut().insert(key.to_string(), value);
    }

    fn get_attrs(&self) -> &HashMap<String, String>;

    fn get_attrs_mut(&mut self) -> &mut HashMap<String, String> {
        unimplemented!()
    }

    fn bed(&self) -> String;
    fn gtf(&self) -> String;
    fn gff(&self) -> String;

    // len and overlaps/contains are now automatically available from EntryT's interval
    fn len(&self) -> usize {
        (self.interval().end - self.interval().start) + 1
    }

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
    pub g_id: Option<String>,
    pub g_parent_id: Option<String>,
    pub attrs: HashMap<String, String>,
    extra_attrs: HashMap<String,String>, // extra attributes that are not part of the GFF/GTF 9th column
    
    pub interval: Interval<usize>,
    children: ArrayBackedIntervalTree<GffObject>, // TODO: add children - any object should be able to inherit children when converted from something like transcript or gene, or anyhting else. It might not be able to do anything with the children, but should be able to store them
    indexed: bool,
    parent: Option<NonNull<GffObject>>, // Pointer to the parent
}

impl Default for GffObject {
    fn default() -> Self {
        GffObject {
            seqid: String::new(),
            source: String::from("GANLIB"),
            g_type: Types::Unknown,
            g_id: None,
            g_parent_id: None,
            strand: '.',
            attrs: HashMap::new(),
            extra_attrs: HashMap::new(),

            interval: Interval::new(0..0).unwrap(),
            children: ArrayBackedIntervalTree::new(),
            indexed: false,
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
            (obj.g_id,obj.g_parent_id) = extract_ids(&obj.attrs, &obj.g_type);

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
    fn new(line: &str) -> Result<GffObject, Box<dyn Error>> {
        GffObject::try_from(line)
    }
    fn parent(&self) -> Option<&dyn GffObjectT> {
        self.parent.and_then(|parent_ptr| unsafe { Some(parent_ptr.as_ref() as &dyn GffObjectT) })
    }
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
    fn get_attrs_mut(&mut self) -> &mut HashMap<String, String> {
        &mut self.attrs
    }
}

impl GffObject {
    pub fn new(line: &str) -> Result<GffObject, Box<dyn Error>> {
        GffObject::try_from(line)
    }

    pub fn add(&mut self, mut obj: GffObject) -> Result<usize, Box<dyn Error>> {
        // Set the child's parent to self
        obj.parent = Some(NonNull::new(self).unwrap());

        self.children.insert(obj);
        self.indexed = false;

        Ok(self.children.len() - 1)
    }

    pub fn get(&self, id: usize) -> Option<&GffObject> {
        self.children.get(id)
    }

    pub fn index(&mut self) {
        if !self.indexed {
            self.children.index();
            self.indexed = true;
        }
    }

    pub fn get_parent(&self) -> Option<&GffObject> {
        self.parent.and_then(|parent_ptr| unsafe { Some(parent_ptr.as_ref()) })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_gffobject_from_default() {
        let mut obj = GffObject::default();
        assert_eq!(obj.seqid, "");
        assert_eq!(obj.source, "GANLIB");
        assert_eq!(obj.interval, Interval::new(0..0).unwrap());
        assert_eq!(obj.strand, '.');
        assert_eq!(obj.attrs.len(), 0);
        assert_eq!(obj.g_type, Types::Unknown);
        assert_eq!(obj.extra_attrs.len(), 0);
        obj.seqid = "chr1".to_string();
        obj.source = "test".to_string();
        obj.interval = Interval::new(1..10).unwrap();
        obj.strand = '+';
        obj.attrs.insert("test".to_string(), "test".to_string());
        obj.g_type = Types::Gene;
        obj.extra_attrs.insert("test".to_string(), "test".to_string());
        assert_eq!(obj.seqid, "chr1");
        assert_eq!(obj.source, "test");
        assert_eq!(obj.interval, Interval::new(1..10).unwrap());
        assert_eq!(obj.strand, '+');
        assert_eq!(obj.attrs.len(), 1);
        assert_eq!(obj.attrs.get("test").unwrap(), "test");
        assert_eq!(obj.g_type, Types::Gene);
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
        assert_eq!(obj.g_type, Types::Gene);
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

    // test iteration over children of the object
    #[test]
    fn test_gffobject_children() {
        let line = "chr1\t.\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\";";
        let mut obj = GffObject::new(line).unwrap();
        let line = "chr1\t.\texon\t11869\t12227\t.\t+\t.\tgene_id \"ENSG00000223972.5\";";
        let exon1 = GffObject::new(line).unwrap();
        let line = "chr1\t.\texon\t12613\t12721\t.\t+\t.\tgene_id \"ENSG00000223972.5\";";
        let exon2 = GffObject::new(line).unwrap();
        let line = "chr1\t.\texon\t13220\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\";";
        let exon3 = GffObject::new(line).unwrap();
        let _ = obj.add(exon1.clone());
        let _ = obj.add(exon2.clone());
        let _ = obj.add(exon3.clone());
        assert_eq!(obj.children.len(), 3);
        let mut iter = obj.children.into_iter();
        assert_eq!(iter.next().unwrap(), &exon1);
        assert_eq!(&exon2, iter.next().unwrap());
        assert_eq!(iter.next().unwrap(), &exon3);

        // iterate and print
        obj.children.into_iter().for_each(|x| println!("{:?}", x.gtf()));
    }

    #[test]
    fn test_types_display() {
        assert_eq!(Types::Gene.to_string(), "gene");
        assert_eq!(Types::Unknown.to_string(), "unknown");
    }

    #[test]
    fn test_default_type() {
        let default = Types::default();
        assert_eq!(default, Types::Unknown);
    }

    #[test]
    fn test_valid_line_parsing() {
        let line = "chr1\tENSEMBL\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=GENE1";
        let obj = GffObject::try_from(line).expect("Failed to parse valid GFF line");
        assert_eq!(obj.seqid, "chr1");
        assert_eq!(obj.interval, Interval::new(100..200).unwrap());
        assert_eq!(obj.strand, '+');
        assert_eq!(obj.g_type, Types::Gene);
        assert_eq!(obj.get_attr("gene_id").unwrap(), "gene1");
    }

    #[test]
    fn test_invalid_line_parsing() {
        let line = "chr1\tENSEMBL\tgene\t100"; // Invalid GFF line
        let obj = GffObject::try_from(line);
        assert!(obj.is_err());
    }

    #[test]
    fn test_obj_search_and_retrieve() {
        // add a few transcripts with exons
        let mut t = GffObject::default();
        let obj = GffObject::new("chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";").unwrap();
        let _ = t.add(obj).unwrap();
        let obj = GffObject::new("chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";").unwrap();
        let _ = t.add(obj).unwrap();
        let obj = GffObject::new("chr1\ttest\texon\t50\t100\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";").unwrap();
        let _ = t.add(obj).unwrap();
        let obj = GffObject::new("chr1\ttest\texon\t1\t40\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";").unwrap();
        let _ = t.add(obj).unwrap();
        t.index();

        // search the interval tree for object
        let objs = t.children.find(Interval::new(1..100).unwrap());
        // print objects
        for obj in objs {
            println!("{:?}", obj);
        }
    }

    #[test]
    fn test_overlaps() {
        let obj1 = GffObject::try_from("chr1\tENSEMBL\tgene\t100\t200\t.\t+\t.\tID=gene1").unwrap();
        let obj2 = GffObject::try_from("chr1\tENSEMBL\texon\t150\t250\t.\t+\t.\tID=exon1").unwrap();
        assert!(obj1.overlaps(&obj2));
        assert!(obj2.overlaps(&obj1));
    }

    #[test]
    fn test_contains() {
        let obj1 = GffObject::try_from("chr1\tENSEMBL\tgene\t100\t200\t.\t+\t.\tID=gene1").unwrap();
        let obj2 = GffObject::try_from("chr1\tENSEMBL\texon\t120\t180\t.\t+\t.\tID=exon1").unwrap();
        assert!(obj1.contains(&obj2));
        assert!(!obj2.contains(&obj1));
    }

    #[test]
    fn test_parent_child_relationship() {
        let mut parent = GffObject::default();
        let child = GffObject::default();
        parent.children.insert(child.clone());
        assert_eq!(parent.children.len(), 1);
    }

    #[test]
    fn test_parent() {
        let mut t = GffObject::new("chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";").unwrap();
        let e1 = GffObject::new("chr1\ttest\texon\t1\t25\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";").unwrap();
        let _ = t.add(e1).unwrap();
        let e2 = GffObject::new("chr1\ttest\texon\t75\t100\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";").unwrap();
        let _ = t.add(e2).unwrap();

        t.index();

        // iterate over the children and check if the parent is correct
        for child in &t.children {
            let parent = child.get_parent().unwrap();
            assert_eq!(parent.interval, t.interval);
        }
    }

    #[test]
    fn test_gff_string_formatting() {
        let line = "chr1\tENSEMBL\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=GENE1";
        let obj = GffObject::try_from(line).unwrap();
        println!("{}", obj.gff());
        let res_line = "chr1\tENSEMBL\tgene\t100\t200\t.\t+\t.\tID=gene1; Name=GENE1;";
        assert_eq!(obj.gff(), res_line);
    }
}
