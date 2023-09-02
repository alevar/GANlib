use std::collections::HashMap;
use bio::data_structures::interval_tree::IntervalTree;
use bio::utils::Interval;
use std::ops::Deref;

use crate::object::*
;

#[derive(Clone, Debug)]
pub struct Exon {
    interval: Interval<u32>,
    source: String,
    tid: Option<String>,
    gid: Option<String>,
    attrs: HashMap<String, String>,
}

impl Exon {
    pub fn new() -> Self {
        Exon {
            interval: Interval::default(),
            source: String::from("GANLIB"),
            tid: None,
            gid: None,
            attrs: HashMap::new(),
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