//! Multi-format sequence file parser.
//!
//! Supports automatic format detection for:
//! - FASTA (.fasta, .fa, .fna, .faa, .fas)
//! - PHYLIP (.phy, .phylip) - sequential and interleaved
//! - NEXUS (.nex, .nexus, .nxs)
//!
//! Format detection priority:
//! 1. Explicit format specification (-f option)
//! 2. File extension
//! 3. Content-based detection

pub mod fasta;
pub mod nexus;
pub mod phylip;

use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use thiserror::Error;

use fasta::parse_fasta_fast;
use crate::model::Alignment;

/// Detected file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileFormat {
    Fasta,
    Phylip,
    Nexus,
}

impl std::fmt::Display for FileFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FileFormat::Fasta => write!(f, "FASTA"),
            FileFormat::Phylip => write!(f, "PHYLIP"),
            FileFormat::Nexus => write!(f, "NEXUS"),
        }
    }
}

/// Errors that can occur during file parsing.
#[derive(Error, Debug)]
pub enum ParseError {
    #[error("Failed to open file: {0}")]
    IoError(#[from] std::io::Error),

    #[error("Empty file")]
    EmptyFile,

    #[error("Could not determine file format.\n\
             Hint: Use -f/--format to specify the format explicitly:\n  \
             seqtui -f fasta <file>   # FASTA format\n  \
             seqtui -f nexus <file>   # NEXUS format\n  \
             seqtui -f phylip <file>  # PHYLIP format")]
    UnknownFormat,

    #[error("Ambiguous file format (could be {possible}).\n\
             Hint: Use -f/--format to specify the format explicitly:\n  \
             seqtui -f {suggestion} <file>")]
    AmbiguousFormat { possible: String, suggestion: String },

    #[error("FASTA error: {0}")]
    FastaError(#[from] fasta::FastaError),

    #[error("PHYLIP error: {0}")]
    PhylipError(#[from] phylip::PhylipError),

    #[error("NEXUS error: {0}")]
    NexusError(#[from] nexus::NexusError),
}

/// Result type for parsing operations.
pub type ParseResult<T> = Result<T, ParseError>;

/// Detects format from file extension.
pub fn detect_format_from_extension<P: AsRef<Path>>(path: P) -> Option<FileFormat> {
    let ext = path.as_ref().extension().and_then(OsStr::to_str)?;
    match ext.to_lowercase().as_str() {
        // FASTA extensions
        "fa" | "fas" | "fasta" | "fna" | "faa" | "ffn" | "frn" => Some(FileFormat::Fasta),
        // NEXUS extensions
        "nex" | "nexus" | "nxs" => Some(FileFormat::Nexus),
        // PHYLIP extensions
        "phy" | "phylip" | "ph" => Some(FileFormat::Phylip),
        _ => None,
    }
}

/// Detects the file format by examining the content.
pub fn detect_format_from_content(content: &str) -> Option<FileFormat> {
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        
        // NEXUS: starts with #NEXUS (case-insensitive) - most specific
        if trimmed.to_uppercase().starts_with("#NEXUS") {
            return Some(FileFormat::Nexus);
        }
        
        // FASTA: starts with > - very clear indicator
        if trimmed.starts_with('>') {
            return Some(FileFormat::Fasta);
        }
        
        // PHYLIP: first line is "ntax nchar" (two integers)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 2 {
            if parts[0].parse::<usize>().is_ok() && parts[1].parse::<usize>().is_ok() {
                return Some(FileFormat::Phylip);
            }
        }
        
        // First non-empty line doesn't match any known format
        return None;
    }
    
    None
}

/// Tries to parse with multiple formats, returning the first success.
fn try_parse_formats(content: &str, formats: &[FileFormat]) -> ParseResult<(Alignment, FileFormat)> {
    let mut last_error = None;
    
    for &format in formats {
        match parse_content(content, format) {
            Ok(alignment) => return Ok((alignment, format)),
            Err(e) => last_error = Some(e),
        }
    }
    
    Err(last_error.unwrap_or(ParseError::UnknownFormat))
}

/// Parses content with a specific format.
fn parse_content(content: &str, format: FileFormat) -> ParseResult<Alignment> {
    match format {
        FileFormat::Fasta => parse_fasta_fast(content).map_err(ParseError::FastaError),
        FileFormat::Phylip => phylip::parse_phylip_str(content).map_err(ParseError::PhylipError),
        FileFormat::Nexus => nexus::parse_nexus_str(content).map_err(ParseError::NexusError),
    }
}

/// Parses a sequence file with optional format specification.
///
/// Detection priority:
/// 1. Explicit format (if provided)
/// 2. File extension
/// 3. Content-based detection
/// 4. Try all formats (FASTA first as most common, then NEXUS, then PHYLIP)
pub fn parse_file_with_options<P: AsRef<Path>>(
    path: P,
    forced_format: Option<FileFormat>,
) -> ParseResult<Alignment> {
    let file = File::open(&path)?;
    let metadata = file.metadata()?;
    let file_size = metadata.len() as usize;
    
    if file_size == 0 {
        return Err(ParseError::EmptyFile);
    }
    
    let mut reader = BufReader::with_capacity(1024 * 1024, file);
    let mut content = String::with_capacity(file_size);
    reader.read_to_string(&mut content)?;
    
    // 1. Use explicit format if provided
    if let Some(format) = forced_format {
        return parse_content(&content, format);
    }
    
    // 2. Try to detect from extension
    if let Some(format) = detect_format_from_extension(&path) {
        match parse_content(&content, format) {
            Ok(alignment) => return Ok(alignment),
            Err(_) => {
                // Extension didn't work, try content detection
            }
        }
    }
    
    // 3. Try content-based detection
    if let Some(format) = detect_format_from_content(&content) {
        match parse_content(&content, format) {
            Ok(alignment) => return Ok(alignment),
            Err(e) => return Err(e),
        }
    }
    
    // 4. Last resort: try all formats in order of likelihood
    // FASTA is most common and has clear markers
    // NEXUS has clear #NEXUS header
    // PHYLIP is most ambiguous
    match try_parse_formats(&content, &[FileFormat::Fasta, FileFormat::Nexus, FileFormat::Phylip]) {
        Ok((alignment, _)) => Ok(alignment),
        Err(_) => Err(ParseError::UnknownFormat),
    }
}

/// Parses a sequence file, automatically detecting the format.
/// Convenience wrapper around parse_file_with_options.
pub fn parse_file<P: AsRef<Path>>(path: P) -> ParseResult<Alignment> {
    parse_file_with_options(path, None)
}

/// Parses a sequence file with explicit format specification.
///
/// Use this when you know the format in advance or want to force a specific parser.
pub fn parse_file_as<P: AsRef<Path>>(path: P, format: FileFormat) -> ParseResult<Alignment> {
    parse_file_with_options(path, Some(format))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_fasta() {
        let content = ">seq1\nACGT\n";
        assert_eq!(detect_format_from_content(content), Some(FileFormat::Fasta));
    }

    #[test]
    fn test_detect_phylip() {
        let content = "  3   10\nseq1      ACGTACGTAC\n";
        assert_eq!(detect_format_from_content(content), Some(FileFormat::Phylip));
    }

    #[test]
    fn test_detect_nexus() {
        let content = "#NEXUS\nBEGIN DATA;\n";
        assert_eq!(detect_format_from_content(content), Some(FileFormat::Nexus));
        
        // Case insensitive
        let content2 = "#nexus\nbegin data;\n";
        assert_eq!(detect_format_from_content(content2), Some(FileFormat::Nexus));
    }

    #[test]
    fn test_detect_unknown() {
        let content = "This is not a valid sequence file\n";
        assert_eq!(detect_format_from_content(content), None);
    }

    #[test]
    fn test_detect_with_leading_empty_lines() {
        let content = "\n\n  \n>seq1\nACGT\n";
        assert_eq!(detect_format_from_content(content), Some(FileFormat::Fasta));
    }

    #[test]
    fn test_detect_from_extension() {
        assert_eq!(detect_format_from_extension("test.fa"), Some(FileFormat::Fasta));
        assert_eq!(detect_format_from_extension("test.fas"), Some(FileFormat::Fasta));
        assert_eq!(detect_format_from_extension("test.fasta"), Some(FileFormat::Fasta));
        assert_eq!(detect_format_from_extension("test.fna"), Some(FileFormat::Fasta));
        assert_eq!(detect_format_from_extension("test.nex"), Some(FileFormat::Nexus));
        assert_eq!(detect_format_from_extension("test.nexus"), Some(FileFormat::Nexus));
        assert_eq!(detect_format_from_extension("test.phy"), Some(FileFormat::Phylip));
        assert_eq!(detect_format_from_extension("test.phylip"), Some(FileFormat::Phylip));
        assert_eq!(detect_format_from_extension("test.txt"), None);
        assert_eq!(detect_format_from_extension("test.aln"), None);
    }

    #[test]
    fn test_parse_real_nexus_file() {
        // Test parsing the actual LOC_01790.nex file if it exists
        let path = "test_data/LOC_01790.nex";
        if std::path::Path::new(path).exists() {
            let result = parse_file(path);
            assert!(result.is_ok(), "Failed to parse LOC_01790.nex: {:?}", result.err());
            let alignment = result.unwrap();
            println!("Found {} sequences", alignment.sequence_count());
            for i in 0..alignment.sequence_count() {
                if let Some(seq) = alignment.get(i) {
                    println!("  {}: {} (len {})", i, seq.id, seq.len());
                }
            }
            assert_eq!(alignment.sequence_count(), 27);
            assert_eq!(alignment.get(0).unwrap().id, "AelongD09");
            assert_eq!(alignment.get(0).unwrap().len(), 860);
        }
    }
}
