use ganlib::prelude::*;


fn main() {
    let mut treader = ganlib::TReader::new(Some("/home/sparrow/genomicTools/ganlib/data/test.gtf")).unwrap();
    treader.add("/home/sparrow/genomicTools/ganlib/data/test2.gtf").unwrap();
    while let Some(gffobj) = treader.next() {
        println!("{:?}", gffobj.gtf());
    }
}