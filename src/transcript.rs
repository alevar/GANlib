use std::collections::HashMap;
use bio::data_structures::interval_tree::IntervalTree;
use bio::utils::Interval;
use std::ops::Deref;

use crate::exon::Exon;

#[derive(Clone, Debug)]
pub struct Transcript {
    source: String,
    exons: IntervalTree<u32,Exon>,
    tid: Option<String>,
    gid: Option<String>,
    attrs: HashMap<String, String>,
}

// impl Deref for Transcript {
//     type Target = Interval;

//     fn deref(&self) -> &Self::Target {
//         &self.interval
//     }
// }

// impl Transcript {
//     pub fn new() -> Self {
//         Transcript {
//             interval: Interval::new(0..0),
//             source: String::from("GANLIB"),
//             obj_type: Types::Transcript,
//             exons: IntervalTree::new(),
//             cds: IntervalTree::new(),
//             tid: None,
//             gid: None,
//             attrs: HashMap::new(),
//             expression: Vec::new(),
//         }
//     }

//     pub fn from_object(obj: Object) -> Self {
//         let mut transcript = obj.to_transcript();
//         transcript.interval = Interval::new(obj.get_start()..obj.get_end());
//         transcript
//     }

//     pub fn clear(&mut self) {
//         self.tid = None;
//         self.gid = None;
//         self.source = String::from("GANLIB");
//         self.obj_type = Types::Transcript;
//         self.interval = Interval::new(0..0);
//         self.exons.clear();
//         self.cds.clear();
//         self.attrs.clear();
//         self.expression.clear();
//     }

//     pub fn add_exon(&mut self, start: u32, end: u32) {
//         self.exons.insert(Interval::new(start..end));
//     }

//     pub fn add_cds(&mut self, start: u32, end: u32) {
//         self.cds.insert(Interval::new(start..end));
//     }

//     pub fn get_exons(&self) -> &IntervalTree {
//         &self.exons
//     }

//     pub fn get_cds(&self) -> &IntervalTree {
//         &self.cds
//     }

//     pub fn get_tid(&self) -> Option<&str> {
//         self.tid.as_deref()
//     }

//     pub fn set_tid(&mut self, tid: String) {
//         self.tid = Some(tid);
//     }

//     pub fn get_gid(&self) -> Option<&str> {
//         self.gid.as_deref()
//     }

//     pub fn set_gid(&mut self, gid: String) {
//         self.gid = Some(gid);
//     }
// }