// defines traits and srtuct for groups of GffObjects

use std::convert::TryFrom;
use std::error::Error;
use std::fs::File;

use bio::utils::Interval;
use bio::data_structures::interval_tree::{ArrayBackedIntervalTree, EntryT};

use std::collections::HashMap;
use std::cmp::Ordering;

use crate::object::{GffObject, GffObjectT};
use crate::transcript::TranscriptRef;
use crate::treader::TReader;
use crate::utils::*;

pub trait GffObjectGroupT {
    type Object: GffObjectT;

    fn new() -> Self;
    fn add_object(&mut self, obj: Self::Object) -> usize;
    fn get(&self, oid: usize) -> Option<&Self::Object>;
    fn get_mut(&mut self, oid: usize) -> Option<&mut Self::Object>;
    fn objects(&self) -> &ArrayBackedIntervalTree<Self::Object>;
    fn objects_mut(&mut self) -> &mut ArrayBackedIntervalTree<Self::Object>;
}

#[derive(Debug)]
pub struct Transcriptome {
    objects: ArrayBackedIntervalTree<GffObject>,
    id_map: HashMap<String, usize>, // map of object IDs to their indices in the tree

    is_indexed: bool,
}

impl GffObjectGroupT for Transcriptome {
    type Object = GffObject;

    fn new() -> Self {
        Transcriptome {
            objects: ArrayBackedIntervalTree::new(),
            id_map: HashMap::new(),
            is_indexed: false,
        }
    }

    fn add_object(&mut self, obj: Self::Object) -> usize {
        self.is_indexed = false;
        // after inserting the object - set its ID to the index in the tree
        // add entry to the id_map if id_str is present
        let oid = self.objects.insert(obj);
        oid
    }

    fn get(&self, oid: usize) -> Option<&Self::Object> {
        Some(self.objects.get(oid).unwrap())
    }

    fn get_mut(&mut self, oid: usize) -> Option<&mut Self::Object> {
        Some(self.objects.get_mut(oid).unwrap())
    }

    fn objects(&self) -> &ArrayBackedIntervalTree<Self::Object> {
        &self.objects
    }
    fn objects_mut(&mut self) -> &mut ArrayBackedIntervalTree<Self::Object> {
        &mut self.objects
    }
}

impl Transcriptome {
    pub fn from_file(fname: &str) -> Result<Transcriptome, Box<dyn Error>> {
        let mut transcriptome = Transcriptome::new();

        let mut reader = TReader::new(Some(fname))?;
        while let Some(obj) = reader.next() {
            transcriptome.add_object(obj);
        }
        transcriptome.is_indexed = false;
        Ok(transcriptome)
    }

    pub fn add_from_file(&mut self, fname: &str) -> Result<(), Box<dyn Error>> {
        let mut reader = TReader::new(Some(fname))?;
        while let Some(obj) = reader.next() {
            self.add_object(obj);
        }
        self.is_indexed = false;
        Ok(())
    }

    fn index(&mut self) {
        // index the tree
        // set is_indexed to true
        self.objects.index();
        self.is_indexed = true;
    }

    fn finalize(&mut self) {
        // finalize the internals of the tree
        // use the attributes of the objects, to figure out the parent/child relationships and set children/parent fields accordingly
        
        for obj in &self.objects {
            let id_str = obj.id_str.clone();
            let parent_id_str = obj.parent_id_str.clone();
            // if parent is not None and parent object exists with the correct ID
            // get parent index and set it as the parent of the current object and add the current object to the parent's children
            match parent_id_str {
                Some(parent_id_str) => {
                    if let Some(parent_oid) = self.id_map.get(&parent_id_str) {
                        let parent_obj = self.objects.get(*parent_oid).unwrap();
                    }
                }
                None => {
                    // if parent is None - and object is of lower types (Exon, CDS, etc)
                    // create a new parent object with the correct ID and add the current object to its children
                    unimplemented!()
                },
            }
        }
        
        
        
        
        
        
        // let attr_id_str = self.objects.get(oid).unwrap().id_str.clone();
        // match attr_id_str {
        //     Some(attr_id_str) => {
        //         self.id_map.insert(attr_id_str, oid);
        //     }
        //     None => (),
        // }
        // // check parent - if present
        // let attr_parent_id_str = self.objects.get(oid).unwrap().parent_id_str.clone();
    }

    pub fn get_transcript<'a>(&'a mut self, tid: usize) -> Option<TranscriptRef<'a, Transcriptome>> {
        Some(TranscriptRef::new(self, tid))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn test_transcriptome() -> Result<(), Box<dyn Error>> {
        let mut transcriptome = Transcriptome::new();
        let obj = GffObject::new("chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";", false)?;
        let oid = transcriptome.add_object(obj);
        Ok(())
    }

    #[test]
    fn test_from_file() {
        let fname = "test.gff";
        let mut file = File::create(fname).unwrap();
        writeln!(file, "##gff-version 3").unwrap();
        writeln!(file, "# comment").unwrap();
        writeln!(file, "chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";").unwrap();
        writeln!(file, "chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";").unwrap();
        file.flush().unwrap();

        let transcriptome = Transcriptome::from_file(fname).unwrap();
        assert_eq!(transcriptome.objects().len(), 2);

        std::fs::remove_file(fname).unwrap();
    }
}