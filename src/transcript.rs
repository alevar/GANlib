use std::collections::HashMap;
use std::fmt::format;
use bio::utils::Interval;
use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use std::cmp::Ordering;
use std::error::Error;

use crate::object::{GffObject, GffObjectT, Types};
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
    extra_attrs: HashMap<String, String>,
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
            extra_attrs: HashMap::new(),
            obj_type: Types::Transcript,
        }
    }
}

impl Transcript{
    fn new(line: &str) -> Result<Self,Box<dyn Error>> {
        match GffObject::try_from(line){
            Ok(obj) => Ok(Transcript::from(obj)),
            Err(e) => Err(e),
        }
    }

    pub fn exons(&self) -> &ArrayBackedIntervalTree<u32,Exon> {
        &self.exons
    }
    pub fn exons_mut(&mut self) -> &mut ArrayBackedIntervalTree<u32,Exon> {
        &mut self.exons
    }

    // get reference to the first exon with a cds entry
    fn first_coding_exon(&self) -> Option<&Exon>{
        match self.exons.into_iter().filter(|x| x.data().is_coding()).next(){
            Some(exon) => Some(exon.data()),
            None => None,
        }
    }
}

impl From<GffObject> for Transcript {
    /// Convert a GffObject into a Transcript
    /// Will only assign transcript_id and gene_id if present in the GffObject
    /// Otherwise if not present will assign None
    /// transcript_id and gene_id need to be assigned manually if not inferred from GffObject
    /// 
    /// # Arguments
    /// 
    /// * `gff_object` - A GffObject
    /// 
    /// # Example
    /// 
    /// ```
    /// use ganlib::prelude::*;
    /// use ganlib::object::GffObject;
    /// use ganlib::transcript::Transcript;
    /// use bio::utils::Interval;
    /// 
    /// let gffobj = GffObject::new("chr1\t.\ttranscript\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972.5\"; transcript_id \"ENST00000456328.2\"; gene_type \"transcribed_unprocessed_pseudogene\"; gene_status \"KNOWN\"; gene_name \"DDX11L1\"; transcript_type \"processed_transcript\"; transcript_status \"KNOWN\"; transcript_name \"DDX11L1-202\"; level 2; havana_gene \"OTTHUMG00000000961.2\"; havana_transcript \"OTTHUMT00000362751.1\"; tag \"basic\"; transcript_support_level \"1\";").unwrap();
    /// let tx = Transcript::from(gffobj);
    /// assert_eq!(tx.seqid(),"chr1");
    /// assert_eq!(tx.strand(),'+');
    /// assert_eq!(tx.interval(),&Interval::new(11869..14409).unwrap());
    /// assert_eq!(tx.source(),".");
    /// assert_eq!(tx.get_attr("gene_id"),Some("ENSG00000223972.5"));
    /// assert_eq!(tx.get_attr("transcript_id"),Some("ENST00000456328.2"));
    /// ```
    fn from(gff_object: GffObject) -> Self {
        let mut tx = Transcript {
                                seqid: gff_object.seqid().to_string(),
                                strand: gff_object.strand(),
                                interval: gff_object.interval().clone(),
                                exons: ArrayBackedIntervalTree::new(),
                                source: gff_object.source().to_string(),
                                obj_type: Types::Exon,
                                attrs: gff_object.get_attrs().clone(),
                                extra_attrs: HashMap::new(),
                                tid: None,
                                gid: None,
                            };

        // extract tid and gid
        tx.tid = tx.get_attr("transcript_id").map(|x| x.to_string());
        // repeat for gid
        tx.gid = tx.get_attr("gene_id").map(|x| x.to_string());
        tx
    }
}

impl GffObjectT for Transcript {
    fn get_type(&self) -> Types {
        Types::Transcript
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
    fn get_attrs(&self) -> &HashMap<String, String> {
        &self.attrs
    }
    fn get_attrs_mut(&mut self) -> &mut HashMap<String, String> {
        &mut self.attrs
    }
}

impl Transcript{
    fn finalize(&mut self){
        self.exons.index();
        self.interval = Interval::new(self.exons.first().interval().start..self.exons.last().interval().end).unwrap();
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