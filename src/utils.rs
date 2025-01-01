use std::collections::HashMap;
use std::fmt::{Formatter, Display};

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Types {
    Gene,
    Transcript,
    Exon,
    CDS,
    UTR,
    Intron,
    Intergenic,
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
            Types::Unknown => "unknown",
        })
    }
}

impl Default for Types {
    fn default() -> Self {
        Types::Unknown
    }
}

pub fn extract_attributes(attr_str: &str) -> HashMap<String, String> {
    let mut attrs = HashMap::new();

    let is_gff = attr_str.contains('=');

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

pub fn extract_ids(attrs: &HashMap<String, String>, feature_type: &Types) -> (Option<String>, Option<String>) {
    // if GFF
    if attrs.contains_key("id") {
        let id = attrs.get("id").cloned();
        let parent = attrs.get("parent").cloned();
        return (id, parent);
    }

    // Else must be GTF
    let id_key = match feature_type {
        Types::Gene => Some("gene_id"),
        Types::Transcript => Some("transcript_id"),
        _ => None,
    };
    let parent_id_key = match feature_type {
        Types::Gene => None,
        Types::Transcript => Some("gene_id"),
        Types::Exon | Types::CDS | Types::Intron | Types::UTR => Some("transcript_id"),
        _ => Some("transcript_id"),
    };

    let id = id_key.and_then(|k| attrs.get(k).cloned());
    let parent_id = parent_id_key.and_then(|k| attrs.get(k).cloned());

    (id, parent_id)
}

fn cleanup_attributes(attrs: &mut HashMap<String, String>) {
    let keys: Vec<String> = attrs.keys().map(|k| k.to_string()).collect();
    for key in keys {
        if key == "gene_id" || key == "transcript_id" || key == "id" || key == "parent" {
            attrs.remove(&key);
        }
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_attributes() {
        let gff_line = "ID=gene1; gene_name=GENE1";
        let gff_attrs = extract_attributes(gff_line);
        assert_eq!(gff_attrs.len(), 2);

        let gtf_line = "gene_id \"gene1\"; gene_name \"GENE1\"";
        let gtf_attrs = extract_attributes(gtf_line);
        assert_eq!(gtf_attrs.len(), 2);
    }

    #[test]
    fn test_extract_ids() {
        let gff_line = "ID=gene1; gene_name=GENE1";
        let gff_attrs = extract_attributes(gff_line);
        let (id, parent) = extract_ids(&gff_attrs, &Types::Gene);
        assert_eq!(id, Some("gene1".to_string()));
        assert_eq!(parent, None);

        let gtf_line = "gene_id \"g1\"; transcript_id \"t1\"";
        let gtf_attrs = extract_attributes(gtf_line);
        let (id, parent) = extract_ids(&gtf_attrs, &Types::Transcript);
        assert_eq!(id, Some("t1".to_string()));
        assert_eq!(parent, Some("g1".to_string()));

        let (id, parent) = extract_ids(&gtf_attrs, &Types::Gene);
        assert_eq!(id, Some("g1".to_string()));
        assert_eq!(parent, None);

        let (id, parent) = extract_ids(&gtf_attrs, &Types::Exon);
        assert_eq!(id, None);
        assert_eq!(parent, Some("t1".to_string()));
    }
}