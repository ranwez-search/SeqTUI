//! SeqTUI - Terminal Alignment Viewer
//!
//! A terminal-based viewer for sequence alignments.
//!
//! ## Usage
//!
//! ```bash
//! seqtui <sequence_file>
//! seqtui -f nexus <sequence_file>  # Force format
//! ```
//!
//! ## Supported Formats
//!
//! - FASTA (.fasta, .fa, .fna, .faa, .fas)
//! - PHYLIP (.phy, .phylip)
//! - NEXUS (.nex, .nexus, .nxs)
//!
//! ## Navigation (Vim-style)
//!
//! - `h/j/k/l`: Move left/down/up/right
//! - `w/b/e`: Word navigation
//! - `:q`: Quit
//! - `:h`: Help

// Use jemalloc for better memory management (returns memory to OS)
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use std::io::{self, Write};
use std::path::PathBuf;

use anyhow::Result;
use clap::{Parser, ValueEnum};

use seqtui::controller::run_app_with_loading;
use seqtui::formats::{parse_file_with_options, FileFormat};
use seqtui::genetic_code::GeneticCodes;
use seqtui::model::{Alignment, Sequence, SequenceType};

/// Runs CLI mode: parse file, optionally translate, and write to output.
fn run_cli_mode(
    file_path: &PathBuf,
    forced_format: Option<FileFormat>,
    output: &str,
    translate: bool,
    genetic_code: u8,
    reading_frame: u8,
) -> Result<()> {
    // Parse the input file
    let alignment = parse_file_with_options(file_path, forced_format)?;

    // Translate if requested
    let output_alignment = if translate {
        if alignment.sequence_type != SequenceType::Nucleotide {
            anyhow::bail!("Cannot translate: input is not a nucleotide sequence");
        }

        let codes = GeneticCodes::new();
        let code = codes.get(genetic_code).ok_or_else(|| {
            anyhow::anyhow!("Unknown genetic code: {}", genetic_code)
        })?;

        // Convert reading frame from 1-6 to 0-2 (frame) and handle reverse complement
        // For now, frames 1-3 are forward (0, 1, 2), frames 4-6 would need reverse complement
        let frame = if reading_frame <= 3 {
            (reading_frame - 1) as usize
        } else {
            // TODO: reverse complement for frames 4-6
            anyhow::bail!("Reverse complement frames (4-6) not yet implemented");
        };

        let translated_seqs: Vec<Sequence> = alignment
            .sequences
            .iter()
            .map(|seq| {
                let aa_data = code.translate_sequence(seq.as_bytes(), frame);
                Sequence::from_bytes(seq.id.clone(), aa_data)
            })
            .collect();

        let mut translated = Alignment::new(translated_seqs);
        translated.sequence_type = SequenceType::AminoAcid;
        translated
    } else {
        alignment
    };

    // Write output
    if output == "-" {
        // Write to stdout
        let stdout = io::stdout();
        let mut handle = stdout.lock();
        for seq in &output_alignment.sequences {
            writeln!(handle, ">{}", seq.id)?;
            writeln!(handle, "{}", seq.as_str())?;
        }
    } else {
        // Write to file
        let mut file = std::fs::File::create(output)?;
        for seq in &output_alignment.sequences {
            writeln!(file, ">{}", seq.id)?;
            writeln!(file, "{}", seq.as_str())?;
        }
        eprintln!(
            "Wrote {} sequences to {}",
            output_alignment.sequence_count(),
            output
        );
    }

    Ok(()
    )
}

/// Extracts the matching key from a sequence ID using the delimiter.
/// If delimiter is None, returns the full ID.
/// If delimiter is Some, returns the first field before the delimiter.
fn extract_key(id: &str, delimiter: Option<&str>) -> String {
    match delimiter {
        Some(delim) => id.split(delim).next().unwrap_or(id).to_string(),
        None => id.to_string(),
    }
}

/// Runs concatenation mode: parse multiple files and concatenate sequences by ID.
fn run_concatenation_mode(
    files: &[PathBuf],
    forced_format: Option<FileFormat>,
    output: &str,
    translate: bool,
    genetic_code: u8,
    reading_frame: u8,
    delimiter: Option<&str>,
    gap_char: Option<char>,  // None = no gap filling, Some(c) = fill with c
    partitions_file: Option<&str>,
    force: bool,
) -> Result<()> {
    use std::collections::{HashMap, HashSet};
    
    let codes = GeneticCodes::new();
    let code = codes.get(genetic_code).ok_or_else(|| {
        anyhow::anyhow!("Unknown genetic code: {}", genetic_code)
    })?;
    let frame = if reading_frame <= 3 {
        (reading_frame - 1) as usize
    } else {
        anyhow::bail!("Reverse complement frames (4-6) not yet implemented");
    };

    // Pass 1: Collect all sequence IDs and validate alignments (if supermatrix)
    let mut all_keys: Vec<String> = Vec::new();
    let mut seen_keys: HashSet<String> = HashSet::new();
    let mut key_file_count: HashMap<String, usize> = HashMap::new(); // count files per key
    let mut file_lengths: Vec<usize> = Vec::new(); // alignment length per file
    
    eprintln!("Pass 1: Scanning {} files...", files.len());
    
    for file_path in files {
        let alignment = parse_file_with_options(file_path, forced_format)?;
        
        // Validate alignment if supermatrix mode (gap filling enabled)
        if gap_char.is_some() && !alignment.is_valid_alignment {
            anyhow::bail!(
                "File {} is not a valid alignment (sequences have different lengths). \
                Supermatrix mode requires aligned sequences.",
                file_path.display()
            );
        }
        
        // Get alignment length (after translation if needed)
        let aln_len = if translate {
            if alignment.sequence_type != SequenceType::Nucleotide {
                anyhow::bail!(
                    "Cannot translate {}: not a nucleotide sequence",
                    file_path.display()
                );
            }
            // Translated length
            (alignment.alignment_length().saturating_sub(frame)) / 3
        } else {
            alignment.alignment_length()
        };
        file_lengths.push(aln_len);
        
        // Collect keys and track which files they appear in
        let mut keys_in_this_file: HashSet<String> = HashSet::new();
        for seq in &alignment.sequences {
            let key = extract_key(&seq.id, delimiter);
            if !seen_keys.contains(&key) {
                seen_keys.insert(key.clone());
                all_keys.push(key.clone());
            }
            keys_in_this_file.insert(key);
        }
        // Increment file count for each key seen in this file
        for key in keys_in_this_file {
            *key_file_count.entry(key).or_insert(0) += 1;
        }
    }
    
    // Sort keys alphabetically for canonical output
    all_keys.sort();
    
    let output_count = all_keys.len();
    // Count IDs that appear in only one file (orphans)
    let orphan_count = key_file_count.values().filter(|&&c| c == 1).count();
    
    eprintln!("Found {} output sequence IDs ({} appear in only one file)", 
              output_count, orphan_count);
    
    // Check for suspicious ID matching: if >30% of output IDs are orphans (only in 1 file)
    let orphan_ratio = orphan_count as f64 / output_count as f64;
    if orphan_ratio > 0.30 && !force {
        // Write log file with all output IDs
        let log_path = "seqtui_ids.log";
        let mut log_file = std::fs::File::create(log_path)?;
        use std::io::Write;
        writeln!(log_file, "# Output sequence IDs from concatenation")?;
        writeln!(log_file, "# {} IDs total, {} appear in only one file ({:.1}%)", 
                 output_count, orphan_count, orphan_ratio * 100.0)?;
        writeln!(log_file, "# IDs marked with * appear in only one file")?;
        writeln!(log_file, "#")?;
        for key in &all_keys {
            let count = key_file_count.get(key).unwrap_or(&0);
            if *count == 1 {
                writeln!(log_file, "{}*", key)?;
            } else {
                writeln!(log_file, "{}", key)?;
            }
        }
        
        anyhow::bail!(
            "Suspicious ID matching: {:.0}% of output IDs appear in only one file ({} / {}).\n\
            This often means sequence names don't match across files.\n\
            - Check if you need -d/--delimiter to extract a common prefix\n\
            - List of output IDs written to: {} (orphans marked with *)\n\
            - Use --force to proceed anyway",
            orphan_ratio * 100.0, orphan_count, output_count, log_path
        );
    }
    
    // Pass 2: Build concatenated sequences
    // seq_data: key -> concatenated sequence bytes
    let mut seq_data: HashMap<String, Vec<u8>> = HashMap::new();
    for key in &all_keys {
        seq_data.insert(key.clone(), Vec::new());
    }
    
    // Track partitions
    let mut partitions: Vec<(String, usize, usize)> = Vec::new(); // (name, start, end)
    let mut current_pos: usize = 1; // 1-based for partition file
    
    eprintln!("Pass 2: Concatenating sequences...");
    
    for (file_idx, file_path) in files.iter().enumerate() {
        let alignment = parse_file_with_options(file_path, forced_format)?;
        let expected_len = file_lengths[file_idx];
        
        // Build map of key -> sequence for this file
        let mut file_seqs: HashMap<String, Vec<u8>> = HashMap::new();
        for seq in &alignment.sequences {
            let key = extract_key(&seq.id, delimiter);
            
            // Translate if needed
            let seq_data = if translate {
                code.translate_sequence(seq.as_bytes(), frame)
            } else {
                seq.as_bytes().to_vec()
            };
            
            // Check for duplicate keys in same file
            if file_seqs.contains_key(&key) {
                eprintln!(
                    "Warning: Duplicate key '{}' in file {} (using first occurrence)",
                    key,
                    file_path.display()
                );
                continue;
            }
            
            file_seqs.insert(key, seq_data);
        }
        
        // Append sequences (or gaps) for each known key
        for key in &all_keys {
            if let Some(seq_bytes) = file_seqs.get(key) {
                seq_data.get_mut(key).unwrap().extend_from_slice(seq_bytes);
            } else if let Some(fill_char) = gap_char {
                // Fill with specified gap character
                let gaps = vec![fill_char as u8; expected_len];
                seq_data.get_mut(key).unwrap().extend_from_slice(&gaps);
            }
            // If no gap_char and sequence is missing, we simply don't extend
            // This means sequences may have different final lengths
        }
        
        // Record partition
        let gene_name = file_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");
        let end_pos = current_pos + expected_len - 1;
        partitions.push((gene_name.to_string(), current_pos, end_pos));
        current_pos = end_pos + 1;
        
        eprintln!(
            "  {} : {} sequences, {} sites",
            file_path.display(),
            file_seqs.len(),
            expected_len
        );
    }
    
    // Validate: if not filling gaps, check that all sequences have the same length
    if gap_char.is_none() {
        let total_len: usize = file_lengths.iter().sum();
        for (key, data) in &seq_data {
            if data.len() != total_len {
                eprintln!(
                    "Warning: Sequence '{}' has length {} (expected {}). \
                    Use -s/--supermatrix to fill missing with gaps.",
                    key,
                    data.len(),
                    total_len
                );
            }
        }
    }
    
    // Write output
    let seq_count = all_keys.len();
    if output == "-" {
        let stdout = io::stdout();
        let mut handle = stdout.lock();
        for key in &all_keys {
            let data = seq_data.get(key).unwrap();
            writeln!(handle, ">{}", key)?;
            // SAFETY: sequence data is valid ASCII
            writeln!(handle, "{}", unsafe { std::str::from_utf8_unchecked(data) })?;
        }
    } else {
        let mut file = std::fs::File::create(output)?;
        for key in &all_keys {
            let data = seq_data.get(key).unwrap();
            writeln!(file, ">{}", key)?;
            writeln!(file, "{}", unsafe { std::str::from_utf8_unchecked(data) })?;
        }
        eprintln!("Wrote {} sequences to {}", seq_count, output);
    }
    
    // Write partitions file if requested
    if let Some(part_file) = partitions_file {
        let mut file = std::fs::File::create(part_file)?;
        for (name, start, end) in &partitions {
            writeln!(file, "{} = {}-{}", name, start, end)?;
        }
        eprintln!("Wrote {} partitions to {}", partitions.len(), part_file);
    }
    
    Ok(())
}

/// File format specification for command line
#[derive(Debug, Clone, Copy, ValueEnum)]
enum FormatArg {
    /// FASTA format
    Fasta,
    /// NEXUS format
    Nexus,
    /// PHYLIP format
    Phylip,
    /// Auto-detect from extension and content
    Auto,
}

impl From<FormatArg> for Option<FileFormat> {
    fn from(arg: FormatArg) -> Self {
        match arg {
            FormatArg::Fasta => Some(FileFormat::Fasta),
            FormatArg::Nexus => Some(FileFormat::Nexus),
            FormatArg::Phylip => Some(FileFormat::Phylip),
            FormatArg::Auto => None,
        }
    }
}

/// SeqTUI - A Vim-style terminal viewer for sequence alignments
///
/// When run without -o/--output, opens an interactive TUI viewer.
/// With -o/--output, runs in CLI mode and writes output to file (or stdout with "-").
/// With multiple input files, concatenates sequences by matching IDs.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Sequence file(s) to process (FASTA, PHYLIP, or NEXUS format)
    /// Multiple files will be concatenated by matching sequence IDs
    #[arg(required = true)]
    files: Vec<PathBuf>,

    /// Force a specific file format (overrides auto-detection)
    #[arg(short = 'f', long = "format", value_enum, default_value = "auto")]
    format: FormatArg,

    /// Output file (enables CLI mode). Use "-" for stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<String>,

    /// Translate nucleotide sequences to amino acids
    #[arg(short = 't', long = "translate")]
    translate: bool,

    /// Genetic code for translation (1-33, default: 1 = Standard)
    #[arg(short = 'g', long = "genetic-code", default_value = "1")]
    genetic_code: u8,

    /// Reading frame for translation (1-3, default: 1)
    #[arg(short = 'r', long = "reading-frame", default_value = "1")]
    reading_frame: u8,

    /// Delimiter for sequence ID matching (use first field before delimiter)
    /// Example: -d "_" matches "Human_gene1" with "Human_gene2" on "Human"
    #[arg(short = 'd', long = "delimiter")]
    delimiter: Option<String>,

    /// Fill missing sequences with gap character (default: '-', or -s '?' -s '.')
    #[arg(short = 's', long = "supermatrix", value_name = "CHAR", default_missing_value = "-", num_args = 0..=1)]
    supermatrix: Option<String>,

    /// Write partition file (gene boundaries for phylogenetic analysis)
    #[arg(short = 'p', long = "partitions")]
    partitions: Option<String>,

    /// Force concatenation even if sequence ID matching looks suspicious
    /// (more than 30% of output IDs appear in only one file)
    #[arg(long = "force")]
    force: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let forced_format: Option<FileFormat> = args.format.into();

    // Validate reading frame (1-3)
    if args.reading_frame < 1 || args.reading_frame > 3 {
        anyhow::bail!("Reading frame must be 1-3 (got {})", args.reading_frame);
    }

    // Validate genetic code (1-33, with some gaps)
    if args.genetic_code < 1 || args.genetic_code > 33 {
        anyhow::bail!("Genetic code must be 1-33 (got {})", args.genetic_code);
    }

    // Parse and validate supermatrix gap character
    let gap_char: Option<char> = match &args.supermatrix {
        None => None,
        Some(s) => {
            if s.len() != 1 {
                anyhow::bail!(
                    "-s/--supermatrix requires a single character (got '{}', {} chars)",
                    s, s.len()
                );
            }
            let c = s.chars().next().unwrap();
            if !c.is_ascii() {
                anyhow::bail!("-s/--supermatrix character must be ASCII (got '{}')", c);
            }
            Some(c)
        }
    };

    // Validate: supermatrix/partitions/delimiter require output mode
    if args.output.is_none() {
        if args.supermatrix.is_some() {
            anyhow::bail!("-s/--supermatrix requires -o/--output");
        }
        if args.partitions.is_some() {
            anyhow::bail!("-p/--partitions requires -o/--output");
        }
        if args.delimiter.is_some() && args.files.len() > 1 {
            anyhow::bail!("-d/--delimiter with multiple files requires -o/--output");
        }
    }

    // Validate: supermatrix/partitions only make sense with multiple files
    if args.files.len() == 1 {
        if args.supermatrix.is_some() {
            anyhow::bail!("-s/--supermatrix requires multiple input files");
        }
        if args.partitions.is_some() {
            anyhow::bail!("-p/--partitions requires multiple input files");
        }
    }

    // Multiple files: concatenation mode (requires -o)
    if args.files.len() > 1 {
        let output = args.output.ok_or_else(|| {
            anyhow::anyhow!("Multiple input files require -o/--output for concatenation")
        })?;
        
        run_concatenation_mode(
            &args.files,
            forced_format,
            &output,
            args.translate,
            args.genetic_code,
            args.reading_frame,
            args.delimiter.as_deref(),
            gap_char,
            args.partitions.as_deref(),
            args.force,
        )?;
    } else {
        // Single file mode
        let file_path = &args.files[0];
        
        if let Some(output) = args.output {
            // CLI mode: output to file/stdout
            run_cli_mode(
                file_path,
                forced_format,
                &output,
                args.translate,
                args.genetic_code,
                args.reading_frame,
            )?;
        } else {
            // TUI mode with optional preset translation settings
            run_app_with_loading(
                file_path.clone(),
                forced_format,
                if args.translate || args.genetic_code != 1 || args.reading_frame != 1 {
                    Some((args.genetic_code, args.reading_frame))
                } else {
                    None
                },
            )?;
        }
    }

    Ok(())
}
