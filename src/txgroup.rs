use crate::object::{ GffObjectT };

// we can define ObjGroupT as a trait to be inherited by Transcriptome, Bundle, Gene, etc
// What if objects can optionally inherit either only GffObjectT or ObjGroupT?

pub trait GffObjectGroupT: GffObjectT{ // subtrait of GffObjectT defines a group of GffObjects which also behaves as a GffObjectT
    fn iter(&self) -> Box<dyn Iterator<Item = &Self>>;
    fn add(&mut self, obj: Self);
    fn remove(&mut self, obj: Self);
    fn num_elements(&self) -> usize;
}

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_group_iteration() {
    //     let grp = TXGroup::new();
    //     grp.iter().for_each(|x| println!("{:?}", x));
    // }
}