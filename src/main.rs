use ganlib::prelude::*;

struct MyStruct {
    data: Vec<i32>,
    index: usize,
}

impl MyStruct {
    fn new(data: Vec<i32>) -> MyStruct {
        MyStruct { data, index: 0 }
    }
}

impl<'a> IntoIterator for &'a MyStruct {
    type Item = &'a i32; // This specifies that the iterator yields references to i32.

    type IntoIter = std::slice::Iter<'a, i32>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.iter()
    }
}

fn main() {
    let my_struct = MyStruct::new(vec![1, 2, 3, 4, 5]);

    for item in &my_struct {
        println!("Item: {}", item);
        break;
    }

    for item in &my_struct {
        println!("Item: {}", item);
    }
}


// fn main() {
//     let mut treader = ganlib::TReader::new(Some("/home/sparrow/genomicTools/ganlib/data/test.gtf")).unwrap();
//     treader.add("/home/sparrow/genomicTools/ganlib/data/test2.gtf").unwrap();
//     while let Some(line) = treader.next() {
//         println!("{:?}", line);
//     }
// }