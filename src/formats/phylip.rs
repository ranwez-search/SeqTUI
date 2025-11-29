//! PHYLIP format parser.
//!
//! Supports both sequential and interleaved PHYLIP formats.
//!
//! ## PHYLIP Format
//!
//! The first line contains the number of sequences and the sequence length:
//! ```text
//!  3 10
//! ```
//!
//! ### Sequential Format
//! Each sequence name (10 chars, padded) followed by all its data:
//! ```text
//!  3 10
//! Seq1      ACGTACGTAC
//! Seq2      TGCATGCATG
//! Seq3      AAAACCCCGG
//! ```
//!
//! ### Interleaved Format
//! Names on first block, then data continues in subsequent blocks:
//! ```text
//!  3 20
//! Seq1      ACGTACGTAC
//! Seq2      TGCATGCATG
//! Seq3      AAAACCCCGG
//!
//! GTGTGTGTGT
//! CACACACACA
//! TTTTTTTTTT
//! ```
//!
//! ## Relaxed Parsing
//!
//! This parser is lenient about:
//! - Name length (not strictly 10 characters)
//! - Whitespace formatting
//! - Case sensitivity

use thiserror::Error;

use crate::model::{Alignment, Sequence};

/// Errors that can occur during PHYLIP parsing.
#[derive(Error, Debug)]
pub enum PhylipError {
    #[error("Empty PHYLIP file")]
    EmptyFile,

    #[error("Invalid header: expected 'ntax nchar' (two integers), got '{0}'")]
    InvalidHeader(String),

    #[error("Invalid sequence count in header: '{0}' is not a valid number")]
    InvalidSequenceCount(String),

    #[error("Invalid sequence length in header: '{0}' is not a valid number")]
    InvalidSequenceLength(String),

    #[error("Expected {expected} sequences but found {found}")]
    SequenceCountMismatch { expected: usize, found: usize },

    #[error("Sequence '{name}' has length {found}, expected {expected}")]
    SequenceLengthMismatch {
        name: String,
        expected: usize,
        found: usize,
    },

    #[error("No sequence data found after header")]
    NoSequenceData,

    #[error("Missing sequence data for '{0}' in interleaved block")]
    MissingInterleavedData(String),

    #[error("Line {line}: {message}")]
    ParseError { line: usize, message: String },
}

/// Result type for PHYLIP operations.
pub type PhylipResult<T> = Result<T, PhylipError>;

/// Parses PHYLIP content from a string.
pub fn parse_phylip_str(content: &str) -> PhylipResult<Alignment> {
    let lines: Vec<&str> = content.lines().collect();
    
    if lines.is_empty() {
        return Err(PhylipError::EmptyFile);
    }

    // Find the header line (first non-empty line)
    let (header_idx, header) = lines
        .iter()
        .enumerate()
        .find(|(_, line)| !line.trim().is_empty())
        .ok_or(PhylipError::EmptyFile)?;

    // Parse header: "ntax nchar"
    let header = header.trim();
    let parts: Vec<&str> = header.split_whitespace().collect();
    
    if parts.len() < 2 {
        return Err(PhylipError::InvalidHeader(header.to_string()));
    }

    let ntax: usize = parts[0]
        .parse()
        .map_err(|_| PhylipError::InvalidSequenceCount(parts[0].to_string()))?;
    
    let nchar: usize = parts[1]
        .parse()
        .map_err(|_| PhylipError::InvalidSequenceLength(parts[1].to_string()))?;

    if ntax == 0 {
        return Err(PhylipError::InvalidSequenceCount("0".to_string()));
    }

    // Collect all lines after header (keep empty lines as block separators)
    let data_lines: Vec<&str> = lines
        .iter()
        .skip(header_idx + 1)
        .copied()
        .collect();

    if data_lines.iter().all(|l| l.trim().is_empty()) {
        return Err(PhylipError::NoSequenceData);
    }

    // Parse - we use a unified approach that handles both sequential and interleaved
    let sequences = parse_phylip_data(&data_lines, ntax, nchar)?;

    Ok(Alignment::new(sequences))
}

/// Splits a line into name and sequence parts.
/// 
/// Handles both strict PHYLIP (10-char names) and relaxed formats.
fn split_name_and_sequence(line: &str) -> (Option<String>, String) {
    let line = line.trim();
    
    if line.is_empty() {
        return (None, String::new());
    }
    
    // Try strict 10-character name first (traditional PHYLIP)
    // In strict format, name is exactly first 10 chars, padded with spaces
    if line.len() >= 10 {
        let potential_name = &line[..10];
        let potential_seq = &line[10..];
        
        let name_trimmed = potential_name.trim();
        let seq_chars: String = potential_seq
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect();
        
        // Use strict format if:
        // 1. Name part has no internal spaces (it's a single word potentially padded)
        // 2. There's actual sequence data after the 10 chars
        // 3. Sequence part is valid
        let name_words: Vec<&str> = name_trimmed.split_whitespace().collect();
        if name_words.len() == 1  // Single word name (with padding)
            && !name_trimmed.is_empty() 
            && !seq_chars.is_empty() 
            && seq_chars.chars().all(|c| is_sequence_char(c)) 
        {
            return (Some(name_trimmed.to_string()), seq_chars);
        }
    }
    
    // Try relaxed format: whitespace-separated name and sequence
    if let Some(space_idx) = line.find(|c: char| c.is_whitespace()) {
        let name = line[..space_idx].trim();
        let rest = &line[space_idx..];
        let seq: String = rest.chars().filter(|c| !c.is_whitespace()).collect();
        
        // Check if this looks like a name + sequence
        if !name.is_empty() && !seq.is_empty() && seq.chars().all(|c| is_sequence_char(c)) {
            return (Some(name.to_string()), seq);
        }
    }
    
    // This line might be just sequence data (continuation line)
    let seq: String = line.chars().filter(|c| !c.is_whitespace()).collect();
    if !seq.is_empty() && seq.chars().all(|c| is_sequence_char(c)) {
        return (None, seq);
    }
    
    // If nothing matches, treat the whole line as a potential name with no sequence
    (Some(line.to_string()), String::new())
}

/// Checks if a character is a valid sequence character.
fn is_sequence_char(c: char) -> bool {
    c.is_ascii_alphabetic() || c == '-' || c == '.' || c == '*' || c == '?'
}

/// Parses PHYLIP data lines. Handles both sequential and interleaved formats.
fn parse_phylip_data(
    lines: &[&str],
    ntax: usize,
    nchar: usize,
) -> PhylipResult<Vec<Sequence>> {
    let mut sequences: Vec<(String, Vec<u8>)> = Vec::with_capacity(ntax);
    let mut in_interleaved_continuation = false;
    let mut interleaved_idx = 0;
    
    for line in lines {
        let trimmed = line.trim();
        
        // Empty line indicates block boundary in interleaved format
        if trimmed.is_empty() {
            if sequences.len() == ntax {
                // All sequences collected, this might be interleaved continuation
                in_interleaved_continuation = true;
                interleaved_idx = 0;
            }
            continue;
        }
        
        let (name, seq) = split_name_and_sequence(trimmed);
        
        if in_interleaved_continuation {
            // Interleaved continuation block (after empty line)
            if interleaved_idx < sequences.len() {
                sequences[interleaved_idx].1.extend(seq.into_bytes());
                interleaved_idx += 1;
            }
        } else if let Some(n) = name {
            // Line with a name - either new sequence or continuation with name (interleaved)
            if sequences.len() < ntax {
                // Still collecting initial sequences
                sequences.push((n, seq.into_bytes()));
            } else {
                // Already have ntax sequences - this might be interleaved with names
                if let Some(pos) = sequences.iter().position(|(sn, _)| sn == &n) {
                    sequences[pos].1.extend(seq.into_bytes());
                }
            }
        } else if !seq.is_empty() {
            // Line without name (continuation of previous sequence in sequential format)
            if !sequences.is_empty() {
                sequences.last_mut().unwrap().1.extend(seq.into_bytes());
            }
        }
        
        // Check if we've filled all sequences to nchar - if so, we're done
        if sequences.len() == ntax && nchar > 0 {
            let all_complete = sequences.iter().all(|(_, data)| data.len() >= nchar);
            if all_complete {
                break;
            }
        }
    }
    
    if sequences.is_empty() {
        return Err(PhylipError::NoSequenceData);
    }

    Ok(sequences
        .into_iter()
        .map(|(name, data)| Sequence::from_bytes(name, data))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sequential_simple() {
        let content = " 3 10
Seq1      ACGTACGTAC
Seq2      TGCATGCATG
Seq3      AAAACCCCGG
";
        let alignment = parse_phylip_str(content).unwrap();
        assert_eq!(alignment.sequence_count(), 3);
        assert_eq!(alignment.get(0).unwrap().id, "Seq1");
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGTACGTAC");
        assert_eq!(alignment.get(1).unwrap().id, "Seq2");
        assert_eq!(alignment.get(2).unwrap().id, "Seq3");
    }

    #[test]
    fn test_parse_sequential_multiline() {
        let content = " 2 20
Seq1      ACGTACGTAC
GGGGGGGGGG
Seq2      TGCATGCATG
CCCCCCCCCC
";
        let alignment = parse_phylip_str(content).unwrap();
        assert_eq!(alignment.sequence_count(), 2);
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGTACGTACGGGGGGGGGG");
        assert_eq!(alignment.get(1).unwrap().as_str(), "TGCATGCATGCCCCCCCCCC");
    }

    #[test]
    fn test_parse_interleaved() {
        let content = " 3 20
Seq1      ACGTACGTAC
Seq2      TGCATGCATG
Seq3      AAAACCCCGG

GGGGGGGGGG
CCCCCCCCCC
TTTTTTTTTT
";
        let alignment = parse_phylip_str(content).unwrap();
        assert_eq!(alignment.sequence_count(), 3);
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGTACGTACGGGGGGGGGG");
        assert_eq!(alignment.get(1).unwrap().as_str(), "TGCATGCATGCCCCCCCCCC");
        assert_eq!(alignment.get(2).unwrap().as_str(), "AAAACCCCGGTTTTTTTTTT");
    }

    #[test]
    fn test_parse_relaxed_names() {
        // Names shorter than 10 chars, using whitespace as delimiter
        let content = "3 10
seq1 ACGTACGTAC
seq2 TGCATGCATG
seq3 AAAACCCCGG
";
        let alignment = parse_phylip_str(content).unwrap();
        assert_eq!(alignment.sequence_count(), 3);
        assert_eq!(alignment.get(0).unwrap().id, "seq1");
    }

    #[test]
    fn test_parse_with_gaps() {
        let content = " 2 10
Seq1      ACGT--GTAC
Seq2      TG--TGCATG
";
        let alignment = parse_phylip_str(content).unwrap();
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGT--GTAC");
        assert_eq!(alignment.get(1).unwrap().as_str(), "TG--TGCATG");
    }

    #[test]
    fn test_empty_file() {
        let content = "";
        assert!(matches!(
            parse_phylip_str(content),
            Err(PhylipError::EmptyFile)
        ));
    }

    #[test]
    fn test_invalid_header() {
        let content = "not a valid header
Seq1 ACGT
";
        // "not" can't be parsed as a number, so we get InvalidSequenceCount
        let result = parse_phylip_str(content);
        assert!(result.is_err());
        
        // Test with single token header (truly invalid format)
        let content2 = "invalid
Seq1 ACGT
";
        assert!(matches!(
            parse_phylip_str(content2),
            Err(PhylipError::InvalidHeader(_))
        ));
    }

    #[test]
    fn test_too_few_sequences() {
        let content = " 3 10
Seq1      ACGTACGTAC
Seq2      TGCATGCATG
";
        // Should handle gracefully - return what we have
        let result = parse_phylip_str(content);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().sequence_count(), 2);
    }
}
