use std::fs::File;
use std::error::Error;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};

use crate::object::{GffObject, GffObjectT};

use crate::utils::*;

// single treader - private struct to parse over a single file
// used in TReader to parse over multiple simultaneously
struct STReader {
    fname: String,
    reader: BufReader<File>,
    comments : Vec<(u32,String)>, // (line number, comment)
    is_gff: Option<bool>,
}

impl STReader{
    pub fn new(fname: &str) -> Result<STReader,Box<dyn Error>>{
        let file = File::open(fname)?;
        let reader = BufReader::new(file);

        // read some lines to determine if gtf or gff
        let mut streader = STReader{fname:fname.to_string(),
                            reader,
                            comments:vec![],
                            is_gff:None};
        
        let res = streader._set_gff();
        match res {
            Ok(_) => (),
            Err(e) => {
                println!("Couldn't determine the file format: {}",e);
                return Err(e);
            }
        }
        Ok(streader)
    }

    fn _set_gff(&mut self) -> Result<(),Box<dyn Error>> {
        for line in self.reader.by_ref().lines() {
            let line = line.unwrap();
            if line.starts_with('#') {
                continue;
            }
            let lcs: Vec<&str> = line.trim().split('\t').collect();
            if lcs.len() != 9 {
                break;
            }
            if !["transcript", "exon", "CDS"].contains(&lcs[2]) {
                break;
            }
            if lcs[2] == "transcript" {
                self.is_gff = Some(is_gff(&line)?);
            }
        }

        match self.is_gff {
            Some(_) => (),
            None => {
                return Err("Couldn't determine the file format. No suitable lines found.".into());
            }
        }
    
        self.reader.rewind().unwrap();

        Ok(())
    }

    pub fn is_gff(&mut self) -> bool {
        // if is_gff hasn't been set - set it
        // return is_gff value
        if self.is_gff.is_none() {
            self._set_gff().unwrap();
        }
        self.is_gff.unwrap()
    }
}

impl Iterator for STReader {
    type Item = String;
    fn next(&mut self) -> Option<Self::Item>{
        let mut line = String::new();
        loop {
            match self.reader.read_line(&mut line) {
                Ok(0) => return None,
                Ok(_) => {
                    if !line.starts_with('#') {
                        return Some(line);
                    }
                    line.clear();
                },
                Err(_) => return None,
            }
        }
    }
}


pub struct TReader {
    fnames: Vec<String>,
    readers: Vec<STReader>,
}

impl Default for TReader {
    fn default() -> Self {
        TReader{fnames:vec![],readers:vec![]}
    }
}

impl TReader {
    pub fn new<A>(args: Option<A>) -> Result<TReader,Box<dyn Error>>
        where A: Into<String> + Copy
    {
        let mut t = TReader::default();
        match args {
            Some(a) => {
                t.add(&a.into())?;
            },
            None => (),
        };
        Ok(t)
    }

    pub fn add(&mut self, fname: &str) -> Result<(),Box<dyn Error>>{
        self.fnames.push(fname.to_string().clone());
        let reader = STReader::new(fname)?;
        self.readers.push(reader);

        Ok(())
    }
}

impl Iterator for TReader {
    type Item = GffObject;
    fn next(&mut self) -> Option<Self::Item>{
        // iterate over readers
        for (i,reader) in self.readers.iter_mut().enumerate() {
            // if reader is not empty, return line
            if let Some(l) = reader.next() {
                
                let robj = match GffObject::new(l.as_str(),reader.is_gff()) {
                    Ok(robj) => robj,
                    Err(e) => {panic!("Error parsing line: {}",e);},
                };
                return Some(robj);
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn test_streader() {
        let fname = "test.gff";
        let mut file = File::create(fname).unwrap();
        writeln!(file,"##gff-version 3").unwrap();
        writeln!(file,"# comment").unwrap();
        writeln!(file,"chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tID=transcript1").unwrap();
        writeln!(file,"chr1\ttest\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=transcript1").unwrap();
        writeln!(file,"chr1\ttest\ttranscript\t1\t100\t.\t+\t.\ttranscript_id=transcript2").unwrap();
        writeln!(file,"chr1\ttest\texon\t1\t100\t.\t+\t.\ttranscript_id=transcript2").unwrap();
        file.flush().unwrap();

        let mut reader = STReader::new(fname).unwrap();
        assert_eq!(reader.is_gff(),true);

        // remove file
        std::fs::remove_file(fname).unwrap();
    }

    #[test]
    fn test_treader() {
        let fname = "test.gff";
        let mut file = File::create(fname).unwrap();
        writeln!(file,"##gff-version 3").unwrap();
        writeln!(file,"# comment").unwrap();
        writeln!(file,"chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tID=transcript1").unwrap();
        writeln!(file,"chr1\ttest\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=transcript1").unwrap();
        writeln!(file,"chr1\ttest\ttranscript\t1\t100\t.\t+\t.\ttranscript_id=transcript2").unwrap();
        writeln!(file,"chr1\ttest\texon\t1\t100\t.\t+\t.\ttranscript_id=transcript2").unwrap();
        file.flush().unwrap();

        let fname = "test2.gff";
        let mut file = File::create(fname).unwrap();
        writeln!(file,"##gff-version 3").unwrap();
        writeln!(file,"# comment").unwrap();
        writeln!(file,"chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tID=transcript1").unwrap();
        writeln!(file,"chr1\ttest\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=transcript1").unwrap();
        writeln!(file,"chr1\ttest\ttranscript\t1\t100\t.\t+\t.\ttranscript_id=transcript2").unwrap();
        writeln!(file,"chr1\ttest\texon\t1\t100\t.\t+\t.\ttranscript_id=transcript2").unwrap();
        file.flush().unwrap();


        // remove files
        std::fs::remove_file("test.gff").unwrap();
        std::fs::remove_file("test2.gff").unwrap();
    }

    #[test]
    fn test_streader_with_invalid_format() {
        let fname = "invalid.gff";
        let mut file = File::create(fname).unwrap();
        writeln!(file, "invalid\tline\twithout\tnine\tcolumns").unwrap();
        writeln!(file, "another invalid line").unwrap();
        file.flush().unwrap();

        let reader_result = STReader::new(fname);
        assert!(reader_result.is_err(), "Invalid file format should return an error");

        std::fs::remove_file(fname).unwrap();
    }

    #[test]
    fn test_treader_iterator_behavior() {
        let fname = "iterator_test.gff";
        let mut file = File::create(fname).unwrap();
        writeln!(file, "##gff-version 3").unwrap();
        writeln!(file, "chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tID=transcript1").unwrap();
        writeln!(file, "chr1\ttest\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=transcript1").unwrap();
        file.flush().unwrap();

        let mut treader = TReader::new(Some(fname)).unwrap();
        let mut count = 0;
        for item in treader {
            println!("{:?}", item);
            count += 1;
        }

        assert_eq!(count, 2, "Should iterate over all lines");

        std::fs::remove_file(fname).unwrap();
    }
}