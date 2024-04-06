use ganlib::prelude::*;
use bio::utils::Interval;

#[test]
fn test_main() {
    let mut treader = ganlib::TReader::new(Some("data/test.gtf")).unwrap();
    treader.add("data/test2.gtf").unwrap();
    while let Some(gffobj) = treader.next() {
        println!("{:?}", gffobj.gtf());
    }
}