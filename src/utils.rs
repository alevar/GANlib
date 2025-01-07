use std::collections::HashMap;
use std::fmt::{Formatter, Display};
use std::error::Error;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Types {
    Gene,
    Transcript,
    Exon,
    CDS,
    UTR,
    Intron,
    Intergenic,
    Other,
    Unknown,
}

impl Display for Types {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", match self {
            Types::Gene => "gene",
            Types::Transcript => "transcript",
            Types::Exon => "exon",
            Types::CDS => "CDS",
            Types::UTR => "UTR",
            Types::Intron => "intron",
            Types::Intergenic => "intergenic",
            Types::Other => "other",
            Types::Unknown => "unknown",
        })
    }
}

impl Default for Types {
    fn default() -> Self {
        Types::Unknown
    }
}

pub fn extract_attributes(attr_str: &str, is_gff: bool) -> HashMap<String, String> {
    let mut attrs = HashMap::new();

    if is_gff {
        for pair in attr_str.split(';').map(str::trim).filter(|s| !s.is_empty()) {
            if let Some((key, value)) = pair.split_once('=') {
                attrs.insert(key.to_lowercase().to_string(), value.to_string());
            }
        }
    } else {
        for pair in attr_str.split(';').map(str::trim).filter(|s| !s.is_empty()) {
            let mut parts = pair.splitn(2, ' ');
            if let (Some(key), Some(value)) = (parts.next(), parts.next()) {
                let value = value.trim_matches('"');
                attrs.insert(key.to_lowercase().to_string(), value.to_string());
            }
        }
    }

    attrs
}

pub fn extract_id(attrs: &HashMap<String, String>, feature_type: &Types, is_gff: bool) -> Option<String> {
    if is_gff {
        if let Some(id) = attrs.get("id") {
            return Some(id.clone());
        }
    } else {
        let id_key = match feature_type {
            Types::Gene => Some("gene_id"),
            Types::Transcript => Some("transcript_id"),
            _ => None,
        };
        if let Some(key) = id_key {
            if let Some(id) = attrs.get(key) {
                return Some(id.clone());
            }
        }
    }
    None
}

pub fn extract_parent_id(attrs: &HashMap<String, String>, feature_type: &Types, is_gff: bool) -> Option<String> {
    if is_gff {
        if let Some(parent) = attrs.get("parent") {
            return Some(parent.clone());
        }
    } else {
        let parent_id_key = match feature_type {
            Types::Gene => None,
            Types::Transcript => Some("gene_id"),
            Types::Exon | Types::CDS | Types::Intron | Types::UTR => Some("transcript_id"),
            _ => Some("transcript_id"),
        };
        if let Some(key) = parent_id_key {
            if let Some(parent) = attrs.get(key) {
                return Some(parent.clone());
            }
        }
    }
    None
}

fn get_attr_value(attrs: &HashMap<String, String>, keys: &[&str]) -> Option<String> {
    let mut value = None;
    for key in keys {
        if let Some(v) = attrs.get(*key) {
            // if value is none - set to value
            // otherwise make sure is the same value
            // otherwise return none
            if value.is_none() {
                value = Some(v.clone());
            } else if value != Some(v.clone()) {
                return None;
            }
        }
    }
    value
}

pub fn attr_is_gff(attrs: &str) -> Result<bool,Box<dyn Error>> {
    // remove any trailing whitespace and last ; if present
    // split by ; and remove any leading or trailing whitespace
    let attrs_vec: Vec<&str> = attrs.trim().trim_end_matches(';').split(';').map(str::trim).collect();

    // if any attribute string begins with ID= or Parent=
    // and every attribute string has "="
    // then is GFF
    let is_gff = attrs_vec.
        iter().
        any(|a| a.starts_with("ID=") || a.starts_with("Parent=")) && attrs_vec.iter().all(|a| a.contains('='));

    if is_gff {
        return Ok(true);
    }
    // if not GFF - need to validate GTF
    // make sure each pair has at least one space followed by \" and ending with \"
    let is_gtf = attrs_vec.iter().all(|a| a.contains(" \"") && a.ends_with("\""));
    if is_gtf {
        return Ok(false);
    }
    else{
        return Err("Cannot determine format".into());
    }
}

// given a gtf/gff line - determine format
pub fn is_gff(line: &str) -> Result<bool,Box<dyn Error>> { // return true if GFF, false if GTF
    // if line starts with # - ERror
    if line.starts_with('#') {
       return Err("Cannot determine format from comment line".into());
    }

    // split line by tab
    let lcs: Vec<&str> = line.trim().split('\t').collect();
    // if not 9 columns - return error
    if lcs.len() != 9 {
        return Err("Invalid number of columns".into());
    }

    attr_is_gff(lcs[8])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_attributes() {
        let gff_line = "ID=gene1; gene_name=GENE1";
        let gff_attrs = extract_attributes(gff_line,true);
        assert_eq!(gff_attrs.len(), 2);

        let gtf_line = "gene_id \"gene1\"; gene_name \"GENE1\"";
        let gtf_attrs = extract_attributes(gtf_line,false);
        assert_eq!(gtf_attrs.len(), 2);
    }

    #[test]
    fn test_extract_ids() {
        let gff_line = "ID=gene1; gene_name=GENE1";
        let gff_attrs = extract_attributes(gff_line,true);
        let id = extract_id(&gff_attrs, &Types::Gene, true);
        let parent = extract_parent_id(&gff_attrs, &Types::Gene, true);
        assert_eq!(id, Some("gene1".to_string()));
        assert_eq!(parent, None);

        let gtf_line = "gene_id \"g1\"; transcript_id \"t1\"";
        let gtf_attrs = extract_attributes(gtf_line,false);
        let id = extract_id(&gtf_attrs, &Types::Transcript, false);
        let parent = extract_parent_id(&gtf_attrs, &Types::Transcript, false);
        assert_eq!(id, Some("t1".to_string()));
        assert_eq!(parent, Some("g1".to_string()));

        let id = extract_id(&gtf_attrs, &Types::Gene, false);
        let parent = extract_parent_id(&gtf_attrs, &Types::Gene, false);
        assert_eq!(id, Some("g1".to_string()));
        assert_eq!(parent, None);

        let id = extract_id(&gtf_attrs, &Types::Exon, false);
        let parent = extract_parent_id(&gtf_attrs, &Types::Exon, false);
        assert_eq!(id, None);
        assert_eq!(parent, Some("t1".to_string()));
    }
}