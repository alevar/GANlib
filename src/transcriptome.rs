

#[derive(Clone,Debug)]
pub struct Transcriptome{

}

// when transcriptome is being finalized the following rules need to apply:
// 1. all transcripts must have a gene_id
// 2. all exons must have a transcript_id
//    - if an exon has None for transcript id - it get's promoted to a single-exon transcript
// 3. all CDS must have a transcript_id