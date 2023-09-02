use ganlib::TReader;

fn main() {
    let mut treader = TReader::new(Some("/home/sparrow/genomicTools/ganlib/data/test.gtf")).unwrap();
    treader.add("/home/sparrow/genomicTools/ganlib/data/test2.gtf").unwrap();
    while let Some(line) = treader.next() {
        println!("{:?}", line);
    }
}