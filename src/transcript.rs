use std::collections::HashMap;

use crate::group::GffObjectGroupT;
use crate::object::GffObjectT;
use crate::utils::*;

use bio::utils::Interval;
use bio::data_structures::interval_tree::EntryT;

pub struct TranscriptRef<'a, Group>
where
    Group: GffObjectGroupT,
{
    parent: &'a mut Group,
    tid: usize,
}

impl<'a, Group> TranscriptRef<'a, Group>
where
Group: GffObjectGroupT,
{
    pub fn change_exon_type(&mut self, gtype: Types) {
        // example of how to apply changes to childred on an object
        let children: Vec<usize> = self
            .parent
            .get(self.tid)
            .unwrap()
            .children()
            .to_vec();

        for child_id in children {
            self
            .parent
            .objects_mut()
            .get_mut(child_id)
            .unwrap()
            .set_type(gtype.clone());
        }
    }
}

impl<'a, Group> EntryT for TranscriptRef<'a, Group>
where
    Group: GffObjectGroupT,
{
    type N = usize;

    fn interval(&self) -> &Interval<Self::N> {
        self.parent.get(self.tid).unwrap().interval()
    }
}

impl<'a, Group> GffObjectT for TranscriptRef<'a, Group>
where
    Group: GffObjectGroupT,
{
    fn seqid(&self) -> &str {
        self.parent.get(self.tid).unwrap().seqid()
    }

    fn strand(&self) -> char {
        self.parent.get(self.tid).unwrap().strand()
    }

    fn get_type(&self) -> Types {
        self.parent.get(self.tid).unwrap().get_type()
    }

    fn source(&self) -> &str {
        self.parent.get(self.tid).unwrap().source()
    }

    fn get_attrs(&self) -> &HashMap<String, String> {
        self.parent.get(self.tid).unwrap().get_attrs()
    }

    fn set_attr(&mut self, key: &str, value: String) {
        self.parent
            .objects_mut()
            .get_mut(self.tid)
            .unwrap()
            .set_attr(key, value);
    }

    fn children(&self) -> &[usize] {
        self.parent.get(self.tid).unwrap().children()
    }

    fn set_type(&mut self, gtype: Types) {
        self.parent.objects_mut().get_mut(self.tid).unwrap().set_type(gtype);
    }

    fn bed(&self) -> String {
        self.parent.get(self.tid).unwrap().bed()
    }

    fn gtf(&self) -> String {
        self.parent.get(self.tid).unwrap().gtf()
    }

    fn gff(&self) -> String {
        self.parent.get(self.tid).unwrap().gff()
    }

    fn score(&self) -> Option<f32> {
        self.parent.get(self.tid).unwrap().score()
    }

    fn phase(&self) -> Option<u32> {
        self.parent.get(self.tid).unwrap().phase()
    }
}

impl<'a, Group> TranscriptRef<'a, Group>
where
    Group: GffObjectGroupT,
{
    pub fn new(parent: &'a mut Group, tid: usize) -> Self {
        TranscriptRef { parent, tid }
    }
}

impl<'a, Group> std::fmt::Display for TranscriptRef<'a, Group>
where
    Group: GffObjectGroupT,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "TranscriptRef: {{ tid: {} }}",
            self.tid
        )
    }
}

impl<'a, Group> std::fmt::Debug for TranscriptRef<'a, Group>
where
    Group: GffObjectGroupT,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "TranscriptRef: {{ tid: {} }}",
            self.tid
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{group::Transcriptome, object::GffObject};

    #[test]
    fn test_transcript_ref() {
        let mut parent = Transcriptome::new();
        let tline = "chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";";
        let tid = parent.add_object(GffObject::new(tline).unwrap());
        let eline = "chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id \"test\"; transcript_id \"test\";";
        let eid = parent.add_object(GffObject::new(eline).unwrap());
        let mut tref = parent.get_transcript(tid).unwrap();
        println!("{:?}", tref);
        tref.change_exon_type(Types::Other);
        println!("{:?}", tref);
    }
}