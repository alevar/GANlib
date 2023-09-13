use ganlib::prelude::*;

// Define the two structs
#[derive(Clone, Debug)]
struct StructA {
    data: i32,
}

#[derive(Clone, Debug)]
struct StructB {
    data: i32,
}

// // Implement the Into trait for converting from StructA to StructB
// impl Into<StructB> for StructA {
//     fn into(self) -> StructB {
//         StructB { data: self.data }
//     }
// }

impl From<StructA> for StructB {
    fn from(other: StructA) -> StructB {
        StructB { data: other.data }
    }
}

// Implement the From trait for converting from StructB to StructA
impl From<StructB> for StructA {
    fn from(other: StructB) -> StructA {
        StructA { data: other.data }
    }
}

fn main() {
    let struct_a = StructA { data: 42 };
    let struct_b: StructB = struct_a.clone().into(); // Convert from StructA to StructB
    let struct_a_again: StructA = StructA::from(struct_b.clone()); // Convert from StructB to StructA

    println!("StructA data: {}", struct_a.data);
    println!("StructB data: {}", struct_b.data);
    println!("StructA data again: {}", struct_a_again.data);
}


// fn main() {
//     let mut treader = ganlib::TReader::new(Some("/home/sparrow/genomicTools/ganlib/data/test.gtf")).unwrap();
//     treader.add("/home/sparrow/genomicTools/ganlib/data/test2.gtf").unwrap();
//     while let Some(line) = treader.next() {
//         println!("{:?}", line);
//     }
// }