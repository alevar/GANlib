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
        let oid = self.objects.insert(obj);
        self.objects.get_mut(oid).unwrap().id = Some(oid);
        // add entry to the id_map if id_str is present
        if let Some(id_str) = self.objects.get(oid).unwrap().id_str.clone() {
            self.id_map.insert(id_str, oid);
        }
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

    pub fn create_parent(&mut self, obj: &GffObject) -> Result<usize, Box<dyn Error>> {
        // create a parent object for the given object
        // return the ID of the parent object
        // if the parent object already exists, return its ID
        // if the parent object does not exist, create it and return its ID
        // if the parent object can not be created, return an error

        // make sure the type is compatible with the type of parent ID extracted
        // for example, exon should have transcript_id (but should also check for the available gene_id as well to be propagated upwards)
        // when creating the parent object - can provide it with the current attribtues, and let it extract form them what is needed
        let mut parent_id: Option<usize> = None;
        match obj.g_type {
            Types::Transcript => {
                // create gene object
            }
            Types::Exon | Types::CDS => {
                // create transcript object
            },
            _ => {
                // create parent object based on the type of the object
            }
        }
        match parent_id {
            Some(pid) => Ok(pid),
            None => Err("Parent object could not be created")?,
        }
    }

    pub fn reset_intervals(&mut self) {
        // objects themselves have no colntrol over the intervals of their children, since children are stored as indices only
        // after children have been added to an object
        // we need to make sure parents have intervals that cover all of their children end-to-end
        
        // for each object, get the min(start), max(end) of its children
        // set the interval of the object to cover all of its children
        // propagate further recursively, so that the parents of the object also have intervals that cover all of their children
    }

    fn finalize(&mut self) -> Result<(), Box<dyn Error>> {
        // finalize the internals of the tree
        // use the attributes of the objects, to figure out the parent/child relationships and set children/parent fields accordingly
        
        // we can not borrow and modify the objects at the same time
        // so we will create a new vector of objects with the correct parent/child relationships
        let mut hierarchy_updates: Vec<(usize, GffObject)> = Vec::new(); // TODO: ideally we do not create copies here, since some internals are heavy (attributes, for example)
        // collect parent-child relationships
        for obj in &self.objects {
            // make sure each object already has an ID assigned (should be handled when creating transcriptome)
            if (obj.id.is_none()){
                Err("Object does not have an ID assigned")?;
            }
            if let Some(parent_id_str) = obj.parent_id_str.clone() {
                // check if parent already exists
                if let Some(parent_oid) = self.id_map.get(&parent_id_str) {
                    hierarchy_updates.push((*parent_oid, obj.clone()));
                }
                {
                    // TODO: parent not found but an ID for it exists
                    // need to create parent object
                    // recursively create parent objects based on the type of the current object, reconstructing the necessary components
                    // unimplemented!()
                }
            }
        }

        // Assigning parent/child relationships to the objects
        for (parent_id, child_obj) in hierarchy_updates {
            if let Some(parent_obj) = self.objects.get_mut(parent_id) {
                parent_obj.add_child(&child_obj);
            } else {
                Err("Parent object not found")?;
            }
        }

        Ok(())
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
    fn test_transcriptome() {
        let mut transcriptome = Transcriptome::new();
        match GffObject::new("chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";", false) {
            Ok(obj) => {
                let oid = transcriptome.add_object(obj);
                let tobj = transcriptome.get(oid);
                let assigned_oid = tobj.unwrap().id.unwrap();
                assert_eq!(assigned_oid, oid);
            },
            Err(e) => {
                eprintln!("Error: {}", e);
                assert!(false);
            }
        }
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

        let mut transcriptome = Transcriptome::from_file(fname).unwrap();
        let res = transcriptome.finalize();
        assert!(res.is_ok());
        assert_eq!(transcriptome.objects().len(), 2);

        std::fs::remove_file(fname).unwrap();
    }
}