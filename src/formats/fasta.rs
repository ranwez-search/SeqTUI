//! FASTA file parser.
//!
//! This module handles reading and parsing FASTA format files.
//! It supports both single-line and multi-line sequences.
//!
//! ## FASTA Format
//!
//! ```text
//! >sequence_identifier optional description
//! ACGTACGTACGT...
//! >another_sequence
//! TGCATGCATGCA...
//! ```

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use thiserror::Error;

use crate::model::{Alignment, Sequence};

/// Errors that can occur during FASTA parsing.
#[derive(Error, Debug)]
pub enum FastaError {
    #[error("Failed to open file: {0}")]
    IoError(#[from] std::io::Error),

    #[error("Empty FASTA file")]
    EmptyFile,

    #[error("Invalid FASTA format: {0}")]
    InvalidFormat(String),

    #[error("Sequence without header at line {0}")]
    SequenceWithoutHeader(usize),
}

/// Result type for FASTA operations.
pub type FastaResult<T> = Result<T, FastaError>;

/// Parses a FASTA file and returns an Alignment.
///
/// # Arguments
///
/// * `path` - Path to the FASTA file
///
/// # Returns
///
/// An `Alignment` containing all sequences from the file.
///
/// # Examples
///
/// ```no_run
/// use seqtui::formats::fasta::parse_fasta_file;
///
/// let alignment = parse_fasta_file("sequences.fasta").unwrap();
/// println!("Loaded {} sequences", alignment.sequence_count());
/// ```
pub fn parse_fasta_file<P: AsRef<Path>>(path: P) -> FastaResult<Alignment> {
    let file = File::open(&path)?;
    let metadata = file.metadata()?;
    let file_size = metadata.len() as usize;
    
    // For large files, read entire file into memory at once (faster than line-by-line)
    if file_size > 1_000_000 {
        // > 1MB: read all at once
        let mut reader = BufReader::with_capacity(1024 * 1024, file); // 1MB buffer
        let mut content = String::with_capacity(file_size);
        reader.read_to_string(&mut content)?;
        parse_fasta_fast(&content)
    } else {
        // Small files: use line-by-line for memory efficiency
        let reader = BufReader::new(file);
        parse_fasta(reader)
    }
}

/// Fast FASTA parser that works on a pre-loaded string.
/// Avoids per-line allocations by working with slices and bytes.
pub fn parse_fasta_fast(content: &str) -> FastaResult<Alignment> {
    // Estimate number of sequences (rough: one per 1KB on average for alignments)
    let estimated_seqs = (content.len() / 1000).max(10);
    let mut sequences = Vec::with_capacity(estimated_seqs);
    
    let mut current_id: Option<&str> = None;
    let mut current_seq: Vec<u8> = Vec::new();
    let mut line_number = 0;
    let mut prev_seq_len: usize = 1000; // Track previous sequence length for better allocation

    for line in content.lines() {
        line_number += 1;
        let line = line.trim();

        // Skip empty lines
        if line.is_empty() {
            continue;
        }

        if let Some(header) = line.strip_prefix('>') {
            // Save previous sequence if exists
            if let Some(id) = current_id.take() {
                if !current_seq.is_empty() {
                    prev_seq_len = current_seq.len(); // Remember length for next allocation
                    current_seq.shrink_to_fit(); // Reclaim excess capacity
                    sequences.push(Sequence::from_bytes(id, std::mem::take(&mut current_seq)));
                }
            }

            // Parse new header - take everything before first space as ID
            let id = header.split_whitespace().next().unwrap_or(header);

            if id.is_empty() {
                return Err(FastaError::InvalidFormat(format!(
                    "Empty sequence identifier at line {}",
                    line_number
                )));
            }

            current_id = Some(id);
            // Use previous sequence length as guide (alignments have uniform length)
            current_seq = Vec::with_capacity(prev_seq_len.max(1000));
        } else {
            // Sequence line
            if current_id.is_none() {
                return Err(FastaError::SequenceWithoutHeader(line_number));
            }

            // Fast append: most FASTA lines don't have internal whitespace
            if line.bytes().all(|b| !b.is_ascii_whitespace()) {
                current_seq.extend_from_slice(line.as_bytes());
            } else {
                // Fallback: filter whitespace (rare case)
                current_seq.extend(line.bytes().filter(|b| !b.is_ascii_whitespace()));
            }
        }
    }

    // Don't forget the last sequence
    if let Some(id) = current_id {
        if !current_seq.is_empty() {
            current_seq.shrink_to_fit(); // Reclaim excess capacity
            sequences.push(Sequence::from_bytes(id, current_seq));
        }
    }

    if sequences.is_empty() {
        return Err(FastaError::EmptyFile);
    }

    sequences.shrink_to_fit(); // Reclaim excess capacity on the vector itself
    Ok(Alignment::new(sequences))
}

/// Parses FASTA content from a reader.
///
/// This function handles both single-line and multi-line sequences.
pub fn parse_fasta<R: BufRead>(reader: R) -> FastaResult<Alignment> {
    let mut sequences = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();
    let mut line_number = 0;
    let mut prev_seq_len: usize = 1000; // Track previous sequence length for better allocation

    for line_result in reader.lines() {
        line_number += 1;
        let line = line_result?;
        let line = line.trim();

        // Skip empty lines
        if line.is_empty() {
            continue;
        }

        if line.starts_with('>') {
            // Save previous sequence if exists
            if let Some(id) = current_id.take() {
                if !current_seq.is_empty() {
                    prev_seq_len = current_seq.len(); // Remember length for next allocation
                    current_seq.shrink_to_fit(); // Reclaim excess capacity
                    sequences.push(Sequence::from_bytes(id, std::mem::take(&mut current_seq)));
                }
            }

            // Parse new header - take everything after '>' and before first space as ID
            let header = &line[1..];
            let id = header
                .split_whitespace()
                .next()
                .unwrap_or(header)
                .to_string();

            if id.is_empty() {
                return Err(FastaError::InvalidFormat(format!(
                    "Empty sequence identifier at line {}",
                    line_number
                )));
            }

            current_id = Some(id);
            // Use previous sequence length as guide (alignments have uniform length)
            current_seq = Vec::with_capacity(prev_seq_len.max(1000));
        } else {
            // Sequence line
            if current_id.is_none() {
                return Err(FastaError::SequenceWithoutHeader(line_number));
            }

            // Append sequence data (removing any whitespace)
            if line.bytes().all(|b| !b.is_ascii_whitespace()) {
                current_seq.extend_from_slice(line.as_bytes());
            } else {
                current_seq.extend(line.bytes().filter(|b| !b.is_ascii_whitespace()));
            }
        }
    }

    // Don't forget the last sequence
    if let Some(id) = current_id {
        if !current_seq.is_empty() {
            current_seq.shrink_to_fit(); // Reclaim excess capacity
            sequences.push(Sequence::from_bytes(id, current_seq));
        }
    }

    if sequences.is_empty() {
        return Err(FastaError::EmptyFile);
    }

    sequences.shrink_to_fit(); // Reclaim excess capacity on the vector itself
    Ok(Alignment::new(sequences))
}

/// Parses FASTA content from a string.
///
/// Useful for testing or processing in-memory data.
pub fn parse_fasta_str(content: &str) -> FastaResult<Alignment> {
    parse_fasta_fast(content)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_fasta() {
        let content = ">seq1\nACGT\n>seq2\nTGCA\n";
        let alignment = parse_fasta_str(content).unwrap();

        assert_eq!(alignment.sequence_count(), 2);
        assert_eq!(alignment.get(0).unwrap().id, "seq1");
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGT");
        assert_eq!(alignment.get(1).unwrap().id, "seq2");
        assert_eq!(alignment.get(1).unwrap().as_str(), "TGCA");
    }

    #[test]
    fn test_parse_multiline_sequence() {
        let content = ">seq1\nACGT\nTGCA\nAAAA\n";
        let alignment = parse_fasta_str(content).unwrap();

        assert_eq!(alignment.sequence_count(), 1);
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGTTGCAAAAA");
    }

    #[test]
    fn test_parse_with_description() {
        let content = ">seq1 This is a description\nACGT\n";
        let alignment = parse_fasta_str(content).unwrap();

        assert_eq!(alignment.get(0).unwrap().id, "seq1");
    }

    #[test]
    fn test_parse_with_empty_lines() {
        let content = ">seq1\nACGT\n\n>seq2\n\nTGCA\n";
        let alignment = parse_fasta_str(content).unwrap();

        assert_eq!(alignment.sequence_count(), 2);
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGT");
        assert_eq!(alignment.get(1).unwrap().as_str(), "TGCA");
    }

    #[test]
    fn test_empty_file() {
        let content = "";
        let result = parse_fasta_str(content);
        assert!(matches!(result, Err(FastaError::EmptyFile)));
    }

    #[test]
    fn test_sequence_without_header() {
        let content = "ACGT\n>seq1\nTGCA\n";
        let result = parse_fasta_str(content);
        assert!(matches!(result, Err(FastaError::SequenceWithoutHeader(_))));
    }

    #[test]
    fn test_alignment_validation() {
        // Valid alignment
        let content = ">seq1\nACGT\n>seq2\nTGCA\n";
        let alignment = parse_fasta_str(content).unwrap();
        assert!(alignment.is_valid_alignment);

        // Invalid alignment (different lengths)
        let content = ">seq1\nACGT\n>seq2\nTG\n";
        let alignment = parse_fasta_str(content).unwrap();
        assert!(!alignment.is_valid_alignment);
        assert!(alignment.warning.is_some());
    }

    #[test]
    fn test_uppercase_preservation() {
        let content = ">seq1\nacgt\n";
        let alignment = parse_fasta_str(content).unwrap();
        // Preserves case as-is
        assert_eq!(alignment.get(0).unwrap().as_str(), "acgt");
    }
}
