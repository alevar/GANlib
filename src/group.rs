// defines traits and srtuct for groups of GffObjects

use std::convert::TryFrom;
use std::error::Error;

use bio::utils::Interval;
use bio::data_structures::interval_tree::{ArrayBackedIntervalTree, EntryT};

use std::collections::HashMap;
use std::cmp::Ordering;

use crate::object::{GffObject, GffObjectT};
use crate::transcript::TranscriptRef;
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
}

impl GffObjectGroupT for Transcriptome {
    type Object = GffObject;

    fn new() -> Self {
        Transcriptome {
            objects: ArrayBackedIntervalTree::new(),
        }
    }

    fn add_object(&mut self, obj: Self::Object) -> usize {
        self.objects.insert(obj);
        self.objects.len() - 1
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
    pub fn get_transcript<'a>(&'a mut self, tid: usize) -> Option<TranscriptRef<'a, Transcriptome>> {
        Some(TranscriptRef::new(self, tid))
    }
}