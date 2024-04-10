use crate::object::{ GffObjectT, GffObject };
use crate::Exon;

// we can define ObjGroupT as a trait to be inherited by Transcriptome, Bundle, Gene, etc
// What if objects can optionally inherit either only GffObjectT or ObjGroupT?
pub trait GffObjectGroupT: GffObjectT{ // subtrait of GffObjectT defines a group of GffObjects which also behaves as a GffObjectT
    type Child: GffObjectT;

    fn add<T>(&mut self, obj: T)
    where
        T: GffObjectT + Into<Exon> + Into<GffObject>;

        fn iter(&self) -> impl Iterator<Item = &Self::Child> + '_;
    fn num_elements(&self) -> usize;
}

#[derive(Clone,Debug)]
pub struct Bundle{

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