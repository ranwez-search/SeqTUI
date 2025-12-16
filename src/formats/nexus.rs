//! NEXUS format parser.
//!
//! Supports the DATA and CHARACTERS blocks commonly used for sequence alignments.
//!
//! ## NEXUS Format
//!
//! NEXUS files start with `#NEXUS` and contain blocks:
//! ```text
//! #NEXUS
//! BEGIN DATA;
//!   DIMENSIONS NTAX=3 NCHAR=10;
//!   FORMAT DATATYPE=DNA GAP=- MISSING=?;
//!   MATRIX
//!     seq1 ACGTACGTAC
//!     seq2 TGCATGCATG
//!     seq3 AAAACCCCGG
//!   ;
//! END;
//! ```
//!
//! ## Supported Features
//!
//! - DATA and CHARACTERS blocks
//! - DIMENSIONS command (NTAX, NCHAR)
//! - FORMAT command (DATATYPE, GAP, MISSING, INTERLEAVE)
//! - MATRIX command (sequential and interleaved)
//!
//! ## Relaxed Parsing
//!
//! - Case insensitive commands
//! - Flexible whitespace
//! - Quoted and unquoted taxon names

use thiserror::Error;

use crate::model::{Alignment, Sequence};

/// Errors that can occur during NEXUS parsing.
#[derive(Error, Debug)]
pub enum NexusError {
    #[error("Not a NEXUS file (must start with #NEXUS)")]
    NotNexus,

    #[error("Empty NEXUS file")]
    EmptyFile,

    #[error("No DATA or CHARACTERS block found")]
    NoDataBlock,

    #[error("Missing DIMENSIONS command in {block} block")]
    MissingDimensions { block: String },

    #[error("Missing MATRIX command in {block} block")]
    MissingMatrix { block: String },

    #[error("Invalid DIMENSIONS: {0}")]
    InvalidDimensions(String),

    #[error("NTAX not specified in DIMENSIONS")]
    MissingNtax,

    #[error("Expected {expected} sequences (NTAX), found {found}")]
    SequenceCountMismatch { expected: usize, found: usize },

    #[error("Sequence '{name}' has length {found}, expected {expected} (NCHAR)")]
    SequenceLengthMismatch {
        name: String,
        expected: usize,
        found: usize,
    },

    #[error("Unterminated MATRIX (missing ';')")]
    UnterminatedMatrix,

    #[error("Duplicate sequence name: '{0}'")]
    DuplicateName(String),

    #[error("Parse error at line {line}: {message}")]
    ParseError { line: usize, message: String },
}

/// Result type for NEXUS operations.
pub type NexusResult<T> = Result<T, NexusError>;

/// Parses NEXUS content from a string.
pub fn parse_nexus_str(content: &str) -> NexusResult<Alignment> {
    let lines: Vec<&str> = content.lines().collect();
    
    if lines.is_empty() {
        return Err(NexusError::EmptyFile);
    }

    // Verify #NEXUS header
    let first_non_empty = lines
        .iter()
        .find(|line| !line.trim().is_empty())
        .ok_or(NexusError::EmptyFile)?;
    
    if !first_non_empty.trim().to_uppercase().starts_with("#NEXUS") {
        return Err(NexusError::NotNexus);
    }

    // Find DATA or CHARACTERS block and extract its content
    let block_content = find_data_block(&lines)?;
    
    // Parse the block
    parse_data_block(&block_content)
}

/// Finds and extracts the content of a DATA or CHARACTERS block.
fn find_data_block(lines: &[&str]) -> NexusResult<String> {
    let mut in_block = false;
    let mut block_lines: Vec<&str> = Vec::new();

    for line in lines {
        let trimmed = line.trim();
        let upper = trimmed.to_uppercase();

        if !in_block {
            // Look for BEGIN DATA; or BEGIN CHARACTERS;
            if upper.starts_with("BEGIN") && (upper.contains("DATA") || upper.contains("CHARACTERS")) {
                in_block = true;
            }
        } else {
            // Check for END;
            if upper.starts_with("END") || upper.starts_with("ENDBLOCK") {
                break;
            }
            block_lines.push(trimmed);
        }
    }

    if block_lines.is_empty() && !in_block {
        return Err(NexusError::NoDataBlock);
    }

    Ok(block_lines.join("\n"))
}

/// Parses the content of a DATA or CHARACTERS block.
fn parse_data_block(content: &str) -> NexusResult<Alignment> {
    // First, normalize the content by joining multi-line commands
    // NEXUS commands end with ';', so we join lines until we see ';'
    let normalized = normalize_nexus_commands(content);
    
    let mut ntax: Option<usize> = None;
    let mut nchar: Option<usize> = None;
    let mut interleave = false;
    let mut matchchar: Option<char> = None;
    let mut matrix_content = String::new();
    let mut in_matrix = false;

    for line in normalized.lines() {
        let trimmed = line.trim();
        let upper = trimmed.to_uppercase();

        // Skip empty lines (but NOT comments - we need them in the matrix)
        if trimmed.is_empty() {
            continue;
        }
        
        // Skip comments outside of matrix
        if !in_matrix && trimmed.starts_with('[') {
            continue;
        }

        if upper.starts_with("DIMENSIONS") {
            // Parse DIMENSIONS NTAX=n NCHAR=n;
            if let Some(n) = extract_param(&upper, "NTAX") {
                ntax = n.parse().ok();
            }
            if let Some(n) = extract_param(&upper, "NCHAR") {
                nchar = n.parse().ok();
            }
        } else if upper.starts_with("FORMAT") {
            // Check for INTERLEAVE in the full FORMAT command
            interleave = upper.contains("INTERLEAVE");
            // Check for MATCHCHAR
            if let Some(mc) = extract_param(&upper, "MATCHCHAR") {
                matchchar = mc.chars().next();
            }
        } else if upper.starts_with("MATRIX") || in_matrix {
            in_matrix = true;
            // This is matrix content - preserve everything including comments
            if upper.starts_with("MATRIX") {
                // Get content after MATRIX keyword
                if let Some(idx) = upper.find("MATRIX") {
                    let after = trimmed[idx + 6..].trim();
                    if !after.is_empty() && after != ";" {
                        matrix_content.push_str(after.trim_end_matches(';'));
                        matrix_content.push(' ');
                    }
                }
            } else {
                let clean = trimmed.trim_end_matches(';');
                if !clean.is_empty() {
                    matrix_content.push_str(clean);
                    matrix_content.push(' ');
                }
            }
            // Check if this line ends the matrix
            if trimmed.ends_with(';') {
                in_matrix = false;
            }
        }
    }

    // Use ntax and nchar if available
    let expected_ntax = ntax.unwrap_or(0);
    let expected_nchar = nchar.unwrap_or(0);

    // Parse the matrix using token-based approach
    let sequences = parse_matrix(&matrix_content, expected_ntax, expected_nchar, interleave, matchchar)?;

    if sequences.is_empty() {
        return Err(NexusError::NoDataBlock);
    }

    Ok(Alignment::new(sequences))
}

/// Normalizes NEXUS commands by joining multi-line commands into single lines.
/// NEXUS commands end with ';', so lines are joined until a ';' is encountered.
/// The MATRIX command is special - we preserve its internal structure.
fn normalize_nexus_commands(content: &str) -> String {
    let mut result = String::new();
    let mut current_command = String::new();
    let mut in_matrix = false;

    for line in content.lines() {
        let trimmed = line.trim();
        
        // Skip empty lines outside of matrix
        if trimmed.is_empty() {
            if in_matrix {
                result.push('\n');
            }
            continue;
        }
        
        // Inside matrix: preserve comments as markers for new sequences
        if in_matrix {
            // Check if matrix ends
            let line_no_comments = remove_nexus_comments(trimmed);
            if line_no_comments.trim().ends_with(';') {
                // Output the line without the trailing semicolon
                let clean = line_no_comments.trim().trim_end_matches(';');
                if !clean.is_empty() {
                    result.push_str(clean);
                    result.push('\n');
                }
                in_matrix = false;
            } else {
                // Keep original line with comments intact for matrix parsing
                result.push_str(trimmed);
                result.push('\n');
            }
            continue;
        }
        
        // Remove comments from non-matrix lines
        let line_no_comments = remove_nexus_comments(trimmed);
        let line_no_comments = line_no_comments.trim();
        
        if line_no_comments.is_empty() {
            continue;
        }

        let upper = line_no_comments.to_uppercase();
        
        if upper.starts_with("MATRIX") {
            // Start of MATRIX - flush any pending command
            if !current_command.is_empty() {
                result.push_str(&current_command);
                result.push('\n');
                current_command.clear();
            }
            in_matrix = true;
            result.push_str(line_no_comments);
            result.push('\n');
        } else {
            // Regular command - accumulate until ';'
            if !current_command.is_empty() {
                current_command.push(' ');
            }
            current_command.push_str(line_no_comments);
            
            if line_no_comments.ends_with(';') {
                result.push_str(&current_command);
                result.push('\n');
                current_command.clear();
            }
        }
    }

    // Flush any remaining command
    if !current_command.is_empty() {
        result.push_str(&current_command);
        result.push('\n');
    }

    result
}

/// Removes NEXUS comments (bracketed text) from a line.
fn remove_nexus_comments(line: &str) -> String {
    let mut result = String::new();
    let mut in_comment = false;
    
    for c in line.chars() {
        if c == '[' {
            in_comment = true;
        } else if c == ']' {
            in_comment = false;
        } else if !in_comment {
            result.push(c);
        }
    }
    
    result
}

/// Extracts a parameter value from a NEXUS command line.
fn extract_param<'a>(line: &'a str, param: &str) -> Option<&'a str> {
    let idx = line.find(param)?;
    let after = &line[idx + param.len()..];
    
    // Skip to '='
    let eq_idx = after.find('=')?;
    let after_eq = &after[eq_idx + 1..];
    
    // Value ends at whitespace, ';', or end
    let end = after_eq
        .find(|c: char| c.is_whitespace() || c == ';')
        .unwrap_or(after_eq.len());
    
    let value = after_eq[..end].trim();
    if value.is_empty() {
        None
    } else {
        Some(value)
    }
}

/// Parses the MATRIX data using a token-based approach.
/// 
/// NEXUS MATRIX format: tokens are separated by whitespace and comments.
/// In sequential format: name1 seq1 name2 seq2 ... (each seq can span multiple lines)
/// We use NTAX to know how many sequences to expect, and NCHAR to know sequence length.
/// If matchchar is specified, it will be replaced with the corresponding character from the first sequence.
fn parse_matrix(content: &str, ntax: usize, nchar: usize, interleave: bool, matchchar: Option<char>) -> NexusResult<Vec<Sequence>> {
    // Tokenize: split by whitespace, remove comments
    let tokens = tokenize_matrix(content);
    
    if tokens.is_empty() {
        return Ok(Vec::new());
    }
    
    // Parse into (name, data) pairs first
    let mut raw_sequences = if interleave && ntax > 0 {
        parse_interleaved_tokens_raw(&tokens, ntax, nchar)?
    } else {
        parse_sequential_tokens_raw(&tokens, ntax, nchar)?
    };
    
    // Apply MATCHCHAR substitution if specified
    if let Some(mc) = matchchar {
        apply_matchchar_raw(&mut raw_sequences, mc);
    }
    
    // Shrink excess capacity from sequence building
    for (_, data) in &mut raw_sequences {
        data.shrink_to_fit();
    }
    
    // Convert to Sequence objects
    let mut result: Vec<Sequence> = raw_sequences
        .into_iter()
        .map(|(name, data)| Sequence::from_bytes(name, data))
        .collect();
    result.shrink_to_fit();
    Ok(result)
}

/// Apply MATCHCHAR substitution on raw sequence data
fn apply_matchchar_raw(sequences: &mut [(String, Vec<u8>)], matchchar: char) {
    if sequences.len() < 2 {
        return;
    }
    
    let matchchar_byte = matchchar as u8;
    let first_seq_data: Vec<u8> = sequences[0].1.clone();
    
    for (_, seq_data) in sequences.iter_mut().skip(1) {
        for (i, byte) in seq_data.iter_mut().enumerate() {
            if *byte == matchchar_byte && i < first_seq_data.len() {
                *byte = first_seq_data[i];
            }
        }
    }
}

/// Tokenize the matrix content: split by whitespace, remove comments.
fn tokenize_matrix(content: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut current_token = String::new();
    let mut in_comment = false;
    let mut in_single_quote = false;
    let mut in_double_quote = false;
    
    for c in content.chars() {
        if in_comment {
            if c == ']' {
                in_comment = false;
            }
            continue;
        }
        
        if c == '[' && !in_single_quote && !in_double_quote {
            // Start of comment - flush current token
            if !current_token.is_empty() {
                tokens.push(std::mem::take(&mut current_token));
            }
            in_comment = true;
            continue;
        }
        
        if c == '\'' && !in_double_quote {
            in_single_quote = !in_single_quote;
            current_token.push(c);
            continue;
        }
        
        if c == '"' && !in_single_quote {
            in_double_quote = !in_double_quote;
            current_token.push(c);
            continue;
        }
        
        if c.is_whitespace() && !in_single_quote && !in_double_quote {
            if !current_token.is_empty() {
                tokens.push(std::mem::take(&mut current_token));
            }
        } else if c != ';' {
            current_token.push(c);
        }
    }
    
    if !current_token.is_empty() {
        tokens.push(current_token);
    }
    
    tokens
}

/// Parse sequential format tokens into raw (name, data) pairs.
/// Pattern: name1 seqdata1 [seqdata1...] name2 seqdata2 [seqdata2...] ...
/// We use NTAX and NCHAR to determine boundaries.
fn parse_sequential_tokens_raw(tokens: &[String], ntax: usize, nchar: usize) -> NexusResult<Vec<(String, Vec<u8>)>> {
    let mut sequences: Vec<(String, Vec<u8>)> = Vec::new();
    let mut i = 0;
    
    while i < tokens.len() && (ntax == 0 || sequences.len() < ntax) {
        // Next token should be a name
        let name = unquote(&tokens[i]);
        i += 1;
        
        // Collect sequence data until we have NCHAR characters (or see what looks like a new name)
        let mut seq_data = Vec::new();
        
        while i < tokens.len() {
            // If we have NCHAR and we've collected enough, stop
            if nchar > 0 && seq_data.len() >= nchar {
                break;
            }
            
            let token = &tokens[i];
            
            // If we don't have NCHAR info, use heuristics to detect next name:
            // - If we have collected some data and this token looks like a name
            //   (contains non-sequence chars like digits mixed with letters, or is quoted)
            if nchar == 0 && !seq_data.is_empty() && looks_like_name(token) {
                break;
            }
            
            // Add this token as sequence data
            seq_data.extend(token.as_bytes());
            i += 1;
        }
        
        sequences.push((name, seq_data));
    }
    
    Ok(sequences)
}

/// Parse interleaved format tokens into raw (name, data) pairs.
/// In interleaved format, names repeat in each block:
/// name1 seq1_part1 name2 seq2_part1 ... nameN seqN_part1
/// name1 seq1_part2 name2 seq2_part2 ... nameN seqN_part2
/// ...
fn parse_interleaved_tokens_raw(tokens: &[String], ntax: usize, nchar: usize) -> NexusResult<Vec<(String, Vec<u8>)>> {
    let mut sequences: Vec<(String, Vec<u8>)> = Vec::with_capacity(ntax);
    let mut name_to_idx: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut i = 0;
    
    while i < tokens.len() {
        // If we've reached the target length for all sequences, we're done
        if nchar > 0 && sequences.len() == ntax && sequences.iter().all(|(_, data)| data.len() >= nchar) {
            break;
        }
        
        let token = unquote(&tokens[i]);
        i += 1;
        
        // Check if this is a known name (from a previous block)
        if let Some(&idx) = name_to_idx.get(&token) {
            // Append sequence data to existing sequence
            if i < tokens.len() {
                sequences[idx].1.extend(tokens[i].as_bytes());
                i += 1;
            }
        } else if sequences.len() < ntax {
            // New sequence name
            let name = token;
            name_to_idx.insert(name.clone(), sequences.len());
            
            // Get sequence data
            if i < tokens.len() {
                let seq_data = tokens[i].as_bytes().to_vec();
                i += 1;
                sequences.push((name, seq_data));
            } else {
                sequences.push((name, Vec::new()));
            }
        } else {
            // We have all ntax sequences, but this token isn't a known name
            // Treat it as sequence data for the next sequence in round-robin
            let seq_idx = (sequences.len() - 1) % ntax;
            sequences[seq_idx].1.extend(token.as_bytes());
        }
    }
    
    Ok(sequences)
}

/// Remove quotes from a string if present
fn unquote(s: &str) -> String {
    let s = s.trim();
    if (s.starts_with('\'') && s.ends_with('\'')) || (s.starts_with('"') && s.ends_with('"')) {
        s[1..s.len()-1].to_string()
    } else {
        s.to_string()
    }
}

/// Heuristic to detect if a token looks like a sequence name rather than sequence data.
/// This is only used when NCHAR is not available.
fn looks_like_name(token: &str) -> bool {
    // Quoted strings are names
    if token.starts_with('\'') || token.starts_with('"') {
        return true;
    }
    
    // If it contains both letters and digits, it's likely a name (like "seq1", "AelongD09")
    let has_letters = token.chars().any(|c| c.is_ascii_alphabetic());
    let has_digits = token.chars().any(|c| c.is_ascii_digit());
    if has_letters && has_digits {
        return true;
    }
    
    // If it contains underscore, likely a name
    if token.contains('_') {
        return true;
    }
    
    // Short tokens that aren't pure sequence chars might be names
    // But we're conservative here - if in doubt, it's sequence data
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_nexus() {
        let content = r#"#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=3 NCHAR=10;
  FORMAT DATATYPE=DNA GAP=- MISSING=?;
  MATRIX
    seq1 ACGTACGTAC
    seq2 TGCATGCATG
    seq3 AAAACCCCGG
  ;
END;
"#;
        let alignment = parse_nexus_str(content).unwrap();
        assert_eq!(alignment.sequence_count(), 3);
        assert_eq!(alignment.get(0).unwrap().id, "seq1");
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGTACGTAC");
    }

    #[test]
    fn test_parse_interleaved_nexus() {
        let content = r#"#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=2 NCHAR=20;
  FORMAT DATATYPE=DNA INTERLEAVE;
  MATRIX
    seq1 ACGTACGTAC
    seq2 TGCATGCATG

    seq1 GGGGGGGGGG
    seq2 CCCCCCCCCC
  ;
END;
"#;
        let alignment = parse_nexus_str(content).unwrap();
        assert_eq!(alignment.sequence_count(), 2);
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGTACGTACGGGGGGGGGG");
        assert_eq!(alignment.get(1).unwrap().as_str(), "TGCATGCATGCCCCCCCCCC");
    }

    #[test]
    fn test_parse_quoted_names() {
        let content = r#"#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=2 NCHAR=10;
  FORMAT DATATYPE=DNA;
  MATRIX
    'seq 1' ACGTACGTAC
    'seq 2' TGCATGCATG
  ;
END;
"#;
        let alignment = parse_nexus_str(content).unwrap();
        assert_eq!(alignment.get(0).unwrap().id, "seq 1");
        assert_eq!(alignment.get(1).unwrap().id, "seq 2");
    }

    #[test]
    fn test_parse_characters_block() {
        let content = r#"#NEXUS
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=10;
  FORMAT DATATYPE=PROTEIN;
  MATRIX
    seq1 ACDEFGHIKL
    seq2 MNPQRSTVWY
  ;
END;
"#;
        // Note: NTAX not specified, should still work
        let alignment = parse_nexus_str(content).unwrap();
        assert_eq!(alignment.sequence_count(), 2);
    }

    #[test]
    fn test_not_nexus() {
        let content = ">seq1\nACGT\n";
        assert!(matches!(parse_nexus_str(content), Err(NexusError::NotNexus)));
    }

    #[test]
    fn test_no_data_block() {
        let content = r#"#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
END;
"#;
        assert!(matches!(
            parse_nexus_str(content),
            Err(NexusError::NoDataBlock)
        ));
    }

    #[test]
    fn test_case_insensitive() {
        let content = r#"#nexus
begin data;
  dimensions ntax=2 nchar=5;
  format datatype=dna;
  matrix
    seq1 ACGTA
    seq2 TGCAT
  ;
end;
"#;
        let alignment = parse_nexus_str(content).unwrap();
        assert_eq!(alignment.sequence_count(), 2);
    }

    #[test]
    fn test_with_gaps() {
        let content = r#"#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=2 NCHAR=10;
  FORMAT DATATYPE=DNA GAP=-;
  MATRIX
    seq1 ACGT--GTAC
    seq2 TG--TGCATG
  ;
END;
"#;
        let alignment = parse_nexus_str(content).unwrap();
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGT--GTAC");
    }

    #[test]
    fn test_multiline_format() {
        // Test FORMAT command split across multiple lines (like Seaview exports)
        let content = r#"#NEXUS
[saved by seaview on Tue Dec 15 15:49:06 2009]
BEGIN DATA;
  DIMENSIONS NTAX=2 NCHAR=10;
  FORMAT DATATYPE=DNA
  GAP=-
  ;
MATRIX
seq1 ACGT--GTAC
seq2 TG--TGCATG
;
END;
"#;
        let alignment = parse_nexus_str(content).unwrap();
        assert_eq!(alignment.sequence_count(), 2);
        assert_eq!(alignment.get(0).unwrap().id, "seq1");
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGT--GTAC");
    }

    #[test]
    fn test_multiline_sequences() {
        // Test sequences split across multiple lines with inline comments
        // Names have digits/underscores to be clearly distinguishable from sequences
        let content = r#"#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=2 NCHAR=20;
  FORMAT DATATYPE=DNA GAP=-;
MATRIX
[1] seq_1
ACGTACGTAC
GGGGGGGGGG
[2] seq_2
TGCATGCATG
CCCCCCCCCC
;
END;
"#;
        let alignment = parse_nexus_str(content).unwrap();
        assert_eq!(alignment.sequence_count(), 2);
        assert_eq!(alignment.get(0).unwrap().id, "seq_1");
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGTACGTACGGGGGGGGGG");
        assert_eq!(alignment.get(1).unwrap().id, "seq_2");
        assert_eq!(alignment.get(1).unwrap().as_str(), "TGCATGCATGCCCCCCCCCC");
    }

    #[test]
    fn test_matchchar() {
        // Test MATCHCHAR substitution - '.' means same as first sequence
        let content = r#"#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=3 NCHAR=10;
  FORMAT DATATYPE=DNA GAP=- MATCHCHAR=.;
  MATRIX
    seq1 ACGTACGTAC
    seq2 ....TG....
    seq3 T.T.T.T.T.
  ;
END;
"#;
        let alignment = parse_nexus_str(content).unwrap();
        assert_eq!(alignment.sequence_count(), 3);
        
        // seq1: reference sequence
        assert_eq!(alignment.get(0).unwrap().as_str(), "ACGTACGTAC");
        // seq2: ....TG.... -> ACGTTGGTAC (. at positions 0-3,6-9 replaced with seq1)
        assert_eq!(alignment.get(1).unwrap().as_str(), "ACGTTGGTAC");
        // seq3: T.T.T.T.T. -> TCTTTCTTTC (. at odd positions replaced with seq1)
        // Position 0: T, 1: C(seq1[1]), 2: T, 3: T(seq1[3]), 4: T, 5: C(seq1[5]), 6: T, 7: T(seq1[7]), 8: T, 9: C(seq1[9])
        assert_eq!(alignment.get(2).unwrap().as_str(), "TCTTTCTTTC");
    }
}
