use std::collections::HashMap;

use crate::group::GffObjectGroupT;
use crate::object::GffObjectT;
use crate::utils::*;

use bio::utils::Interval;
use bio::data_structures::interval_tree::EntryT;

#[derive(Debug)]
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
            .unwrap();
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::group::Transcriptome;

    #[test]
    fn test_transcript_ref() {
        let mut parent = Transcriptome::new();
    }
}