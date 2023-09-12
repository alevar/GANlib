use std::arch::x86_64::_mm_sha1nexte_epu32;
use std::fs::File;
use std::error::Error;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::str::FromStr;
use noodles_gtf as gtf;
use noodles_gff as gff;
use noodles_gff::Line;

use crate::object::GffObject;
use crate::factory::GTFObjectFactory;

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
        let mut reader = BufReader::new(file);

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
        let mut gff_result = None;

        for line in self.reader.by_ref().lines() {

            let line = line.unwrap();
            if line.starts_with('#') {
                continue;
            }
            let lcs: Vec<&str> = line.trim().split('\t').collect();
            if lcs.len() != 9 {
                gff_result = None;
                break;
            }
            if !["transcript", "exon", "CDS"].contains(&lcs[2]) {
                gff_result = None;
                break;
            }
            if lcs[2] == "transcript" {
                if lcs[8].starts_with("ID=") {
                    gff_result = Some(true);
                    break;
                } else if lcs[8].starts_with("transcript_id") {
                    gff_result = Some(false);
                    break;
                } else {
                    gff_result = None;
                    break;
                }
            }
        }
    
        self.reader.rewind().unwrap();
        
        match gff_result{
            None => Err("Unable to determine file format".into()),
            Some(gff) => {
                self.is_gff = Some(gff);
                Ok(())
            }
        }
    }

    pub fn is_gff(&mut self) -> bool {
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
    current_txs : Vec<Option<String>>,
}

impl Default for TReader {
    fn default() -> Self {
        TReader{fnames:vec![],readers:vec![],current_txs:vec![]}
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
        self.current_txs.push(self.readers.last_mut().unwrap().next());

        Ok(())
    }
}

impl Iterator for TReader {
    type Item = GffObject;
    fn next(&mut self) -> Option<Self::Item>{
        // iterate over readers
        for (i,reader) in self.readers.iter_mut().enumerate() {
            // if reader is empty, skip
            if self.current_txs[i].is_none() {
                continue;
            }
            // if reader is not empty, return line
            if let Some(l) = reader.next() {
                self.current_txs[i] = Some(l.clone());
                let mut gffobj = GffObject::from(l.as_str());
                return Some(gffobj);
            }
            // if reader is empty, set to None
            self.current_txs[i] = None;
        }
        None
    }
}
