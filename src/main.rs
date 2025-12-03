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

use seqtui::controller::{run_app_with_loading, run_app_with_file_browser};
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
    force: bool,
) -> Result<()> {
    // Validate nucleotide content if translating
    if translate {
        validate_nucleotide_files(&[file_path.clone()], forced_format, force, "Translation")?;
    }
    
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

/// Extracts the file basename (without extension) for use as CHROM name
fn get_chrom_name(file_path: &PathBuf) -> String {
    file_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string()
}

/// Computes nucleotide statistics for an alignment.
/// Returns (nt_count, total_count, nt_ratio) where nt_count is the number of
/// ACGT characters and total_count excludes gaps and missing data (N, ?).
fn compute_nt_stats(alignment: &Alignment) -> (usize, usize, f64) {
    let mut nt_count = 0usize;
    let mut total_count = 0usize;
    
    for seq in &alignment.sequences {
        for &b in seq.as_bytes() {
            match b.to_ascii_uppercase() {
                b'A' | b'C' | b'G' | b'T' => {
                    nt_count += 1;
                    total_count += 1;
                }
                b'-' | b'N' | b'?' => {} // Skip gaps and missing
                _ => total_count += 1,   // Non-NT character
            }
        }
    }
    
    let ratio = if total_count > 0 {
        nt_count as f64 / total_count as f64
    } else {
        1.0 // Empty sequences are considered NT
    };
    
    (nt_count, total_count, ratio)
}

/// Validates that files appear to contain nucleotide sequences.
/// Returns Ok(()) if all files pass, or an error with details if not.
/// Files with <50% ACGT characters (excluding gaps/N/?) are flagged.
fn validate_nucleotide_files(
    files: &[PathBuf],
    forced_format: Option<FileFormat>,
    force: bool,
    mode_name: &str,
) -> Result<()> {
    let mut suspect_files: Vec<(String, usize, usize, f64)> = Vec::new();
    
    for file_path in files {
        let alignment = parse_file_with_options(file_path, forced_format)?;
        let (nt_count, total_count, ratio) = compute_nt_stats(&alignment);
        
        if ratio < 0.5 {
            suspect_files.push((
                file_path.display().to_string(),
                nt_count,
                total_count,
                ratio,
            ));
        }
    }
    
    if !suspect_files.is_empty() && !force {
        // Write log file with details
        let log_path = "seqtui_nt_check.log";
        let mut log_file = std::fs::File::create(log_path)?;
        writeln!(log_file, "# Nucleotide content check for {} mode", mode_name)?;
        writeln!(log_file, "# Files with <50% ACGT characters (likely amino acid sequences)")?;
        writeln!(log_file, "#")?;
        writeln!(log_file, "# File\tACGT_count\tTotal_chars\tACGT_ratio")?;
        for (path, nt, total, ratio) in &suspect_files {
            writeln!(log_file, "{}\t{}\t{}\t{:.1}%", path, nt, total, ratio * 100.0)?;
        }
        
        anyhow::bail!(
            "{} mode requires nucleotide sequences, but {} file(s) appear to be amino acids:\n\
            - Less than 50% of characters are ACGT (excluding gaps/N/?)\n\
            - Details written to: {}\n\
            - Use --force to proceed anyway",
            mode_name,
            suspect_files.len(),
            log_path
        );
    }
    
    Ok(())
}

/// Runs VCF mode: extract biallelic SNPs with minimum flanking distance
fn run_vcf_mode(
    files: &[PathBuf],
    forced_format: Option<FileFormat>,
    output: &str,
    min_dist: usize,
    delimiter: Option<&str>,
    force: bool,
) -> Result<()> {
    use std::collections::{HashMap, HashSet};
    
    // Pass 1: Collect all sequence IDs (sorted alphabetically) and validate nucleotide content
    let mut all_keys: Vec<String> = Vec::new();
    let mut seen_keys: HashSet<String> = HashSet::new();
    let mut suspect_files: Vec<(String, usize, usize, f64)> = Vec::new();
    
    eprintln!("Pass 1: Scanning {} file(s) for sequence IDs...", files.len());
    
    for file_path in files {
        let alignment = parse_file_with_options(file_path, forced_format)?;
        
        // Validate: must be a valid alignment
        if !alignment.is_valid_alignment {
            anyhow::bail!(
                "File {} is not a valid alignment (sequences have different lengths). \
                VCF mode requires aligned sequences.",
                file_path.display()
            );
        }
        
        // Validate: must be nucleotide
        if alignment.sequence_type != SequenceType::Nucleotide {
            anyhow::bail!(
                "File {} is not a nucleotide alignment. VCF mode requires nucleotide sequences.",
                file_path.display()
            );
        }
        
        // Check nucleotide content (integrated validation)
        let (nt_count, total_count, ratio) = compute_nt_stats(&alignment);
        if ratio < 0.5 {
            suspect_files.push((
                file_path.display().to_string(),
                nt_count,
                total_count,
                ratio,
            ));
        }
        
        // Collect all keys
        for seq in &alignment.sequences {
            let key = extract_key(&seq.id, delimiter);
            if !seen_keys.contains(&key) {
                seen_keys.insert(key.clone());
                all_keys.push(key);
            }
        }
    }
    
    // Report suspect files if any (unless --force)
    if !suspect_files.is_empty() && !force {
        let log_path = "seqtui_nt_check.log";
        let mut log_file = std::fs::File::create(log_path)?;
        writeln!(log_file, "# Nucleotide content check for VCF mode")?;
        writeln!(log_file, "# Files with <50% ACGT characters (likely amino acid sequences)")?;
        writeln!(log_file, "#")?;
        writeln!(log_file, "# File\tACGT_count\tTotal_chars\tACGT_ratio")?;
        for (path, nt, total, ratio) in &suspect_files {
            writeln!(log_file, "{}\t{}\t{}\t{:.1}%", path, nt, total, ratio * 100.0)?;
        }
        
        anyhow::bail!(
            "VCF mode requires nucleotide sequences, but {} file(s) appear to be amino acids:\n\
            - Less than 50% of characters are ACGT (excluding gaps/N/?)\n\
            - Details written to: {}\n\
            - Use --force to proceed anyway",
            suspect_files.len(),
            log_path
        );
    }
    
    if all_keys.is_empty() {
        anyhow::bail!("No sequences found in input files");
    }
    
    // Sort all keys alphabetically
    all_keys.sort();
    
    eprintln!("Found {} sequences", all_keys.len());
    
    // Prepare VCF output
    let mut vcf_lines: Vec<String> = Vec::new();
    
    // Pass 2: Process each file to find SNPs
    eprintln!("Pass 2: Scanning for biallelic SNPs (min flanking distance: {})...", min_dist);
    
    for file_path in files {
        let alignment = parse_file_with_options(file_path, forced_format)?;
        let chrom = get_chrom_name(file_path);
        let aln_len = alignment.alignment_length();
        
        // Build map: key -> sequence bytes
        let mut seq_map: HashMap<String, &[u8]> = HashMap::new();
        for seq in &alignment.sequences {
            let key = extract_key(&seq.id, delimiter);
            seq_map.insert(key, seq.as_bytes());
        }
        
        // Single pass: collect alleles using bit flags and track if site has only real NTs
        // Bit flags: A=1, C=2, G=4, T=8
        let mut real_nt_only: Vec<bool> = vec![true; aln_len];
        let mut seen_nt: Vec<u8> = vec![0; aln_len];
        
        for seq in &alignment.sequences {
            let seq_bytes = seq.as_bytes();
            for pos in 0..aln_len {
                if real_nt_only[pos] {
                    match seq_bytes[pos].to_ascii_uppercase() {
                        b'A' => seen_nt[pos] |= 1,
                        b'C' => seen_nt[pos] |= 2,
                        b'G' => seen_nt[pos] |= 4,
                        b'T' => seen_nt[pos] |= 8,
                        b'N' | b'?' => {} // Missing data, ignore
                        _ => real_nt_only[pos] = false, // Gap or other -> exclude
                    }
                }
            }
        }
        
        // Derive polymorphic sites: !real_nt_only OR seen_nt has >1 bit set
        // Polymorphic sites reset the flanking distance counter
        let mut reset: Vec<bool> = vec![false; aln_len];
        for pos in 0..aln_len {
            reset[pos] = !real_nt_only[pos] || seen_nt[pos].count_ones() > 1;
        }
        
        // Compute distLeft: distance to nearest polymorphic site on the left
        let mut dist_left: Vec<usize> = vec![0; aln_len];
        for i in 1..aln_len {
            dist_left[i] = if reset[i - 1] { 0 } else { dist_left[i - 1] + 1 };
        }
        
        // Compute distRight: distance to nearest polymorphic site on the right
        let mut dist_right: Vec<usize> = vec![0; aln_len];
        for i in (0..aln_len - 1).rev() {
            dist_right[i] = if reset[i + 1] { 0 } else { dist_right[i + 1] + 1 };
        }
        
        // Write VCF lines for biallelic sites with sufficient flanking distance
        // Biallelic = real_nt_only AND seen_nt has exactly 2 bits set
        let mut snp_count = 0;
        for pos in 0..aln_len {
            // Check biallelic: real_nt_only AND exactly 2 alleles AND sufficient flanking distance
            if !real_nt_only[pos] || seen_nt[pos].count_ones() != 2 
               || dist_left[pos] < min_dist || dist_right[pos] < min_dist {
                continue;
            }
            
            // Decode alleles from bit flags (alphabetical order: A, C, G, T)
            let allele_bits = seen_nt[pos];
            let alleles: Vec<u8> = [b'A', b'C', b'G', b'T']
                .iter()
                .zip([1u8, 2, 4, 8])
                .filter(|(_, bit)| allele_bits & bit != 0)
                .map(|(base, _)| *base)
                .collect();
            
            // REF = first sample (alphabetically) with an actual ACGT at this position
            // ALT = the other allele
            let mut ref_base: Option<u8> = None;
            for key in &all_keys {
                if let Some(seq) = seq_map.get(key) {
                    let b = seq[pos].to_ascii_uppercase();
                    if b == alleles[0] || b == alleles[1] {
                        ref_base = Some(b);
                        break;
                    }
                }
            }
            let ref_base = ref_base.unwrap(); // Safe: biallelic means at least one sample has each allele
            let alt_base = if alleles[0] == ref_base { alleles[1] } else { alleles[0] };
            
            // Build genotypes for all samples
            let gt_strings: Vec<String> = all_keys.iter()
                .map(|key| {
                    if let Some(seq) = seq_map.get(key) {
                        let b = seq[pos].to_ascii_uppercase();
                        if b == ref_base { "0".to_string() }
                        else if b == alt_base { "1".to_string() }
                        else { ".".to_string() } // N, ?, or other -> missing
                    } else {
                        ".".to_string() } // Sequence not in file -> missing
                })
                .collect();
            
            // VCF line: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT samples...
            let vcf_line = format!(
                "{}\t{}\t.\t{}\t{}\t.\tPASS\tDL={};DR={}\tGT\t{}",
                chrom,
                pos + 1, // 1-based position
                ref_base as char,
                alt_base as char,
                dist_left[pos],
                dist_right[pos],
                gt_strings.join("\t")
            );
            vcf_lines.push(vcf_line);
            snp_count += 1;
        }
        
        eprintln!("  {} : {} sites, {} isolated biallelic SNPs selected", chrom, aln_len, snp_count);
    }
    
    // Write VCF output
    let vcf_header = format!(
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=DL,Number=1,Type=Integer,Description=\"Distance to nearest polymorphic site on the left\">\n\
         ##INFO=<ID=DR,Number=1,Type=Integer,Description=\"Distance to nearest polymorphic site on the right\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        all_keys.join("\t")
    );
    
    if output == "-" {
        let stdout = io::stdout();
        let mut handle = stdout.lock();
        writeln!(handle, "{}", vcf_header)?;
        for line in &vcf_lines {
            writeln!(handle, "{}", line)?;
        }
    } else {
        let mut file = std::fs::File::create(output)?;
        writeln!(file, "{}", vcf_header)?;
        for line in &vcf_lines {
            writeln!(file, "{}", line)?;
        }
        eprintln!("Wrote {} SNPs to {}", vcf_lines.len(), output);
    }
    
    Ok(())
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
    
    // Validate nucleotide content if translating
    if translate {
        validate_nucleotide_files(files, forced_format, force, "Translation")?;
    }
    
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
#[command(
    author = "V. Ranwez",
    version,
    about,
    long_about = None,
    after_help = "Repository: https://github.com/ranwez-search/SeqTUI"
)]
struct Args {
    /// Sequence file(s) to process (FASTA, PHYLIP, or NEXUS format)
    /// Multiple files will be concatenated by matching sequence IDs
    /// If omitted, opens file browser in TUI mode
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

    /// Extract biallelic SNPs to VCF with minimum flanking monomorphic distance
    /// Example: -v 300 requires 300 consecutive monomorphic sites on each side
    #[arg(short = 'v', long = "vcf", value_name = "MIN_DIST")]
    vcf: Option<usize>,

    /// Force operation despite warnings (suspicious ID matching or non-nucleotide sequences)
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
        if args.vcf.is_some() {
            anyhow::bail!("-v/--vcf requires -o/--output");
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

    // Validate: VCF mode is incompatible with other CLI options
    if args.vcf.is_some() {
        if args.translate {
            anyhow::bail!("-v/--vcf is incompatible with -t/--translate");
        }
        if args.supermatrix.is_some() {
            anyhow::bail!("-v/--vcf is incompatible with -s/--supermatrix");
        }
        if args.partitions.is_some() {
            anyhow::bail!("-v/--vcf is incompatible with -p/--partitions");
        }
    }

    // No files provided: open TUI with file browser
    if args.files.is_empty() {
        if args.output.is_some() {
            anyhow::bail!("CLI mode (-o/--output) requires at least one input file");
        }
        return run_app_with_file_browser();
    }

    // VCF mode: extract biallelic SNPs
    if let Some(min_dist) = args.vcf {
        let output = args.output.as_ref().unwrap(); // Already validated above
        return run_vcf_mode(
            &args.files,
            forced_format,
            output,
            min_dist,
            args.delimiter.as_deref(),
            args.force,
        );
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
                args.force,
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufRead;
    use std::sync::atomic::{AtomicUsize, Ordering};
    
    static TEST_COUNTER: AtomicUsize = AtomicUsize::new(0);
    
    /// Helper to run VCF mode and capture output lines (excluding header)
    fn run_vcf_test(files: &[&str], min_dist: usize) -> Vec<String> {
        let file_paths: Vec<PathBuf> = files.iter()
            .map(|f| PathBuf::from(format!("test_data/vcf_tests/{}", f)))
            .collect();
        
        // Use unique temp file per test to avoid race conditions
        let test_id = TEST_COUNTER.fetch_add(1, Ordering::SeqCst);
        let tmp_output = format!("/tmp/seqtui_test_vcf_{}.vcf", test_id);
        run_vcf_mode(&file_paths, None, &tmp_output, min_dist, None, false).unwrap();
        
        // Read and return data lines (skip header)
        let file = std::fs::File::open(&tmp_output).unwrap();
        let reader = std::io::BufReader::new(file);
        let lines: Vec<String> = reader.lines()
            .filter_map(|l| l.ok())
            .filter(|l| !l.starts_with('#'))
            .collect();
        
        // Clean up temp file
        let _ = std::fs::remove_file(&tmp_output);
        
        lines
    }
    
    /// Helper to parse a VCF data line
    fn parse_vcf_line(line: &str) -> (String, usize, char, char, usize, usize, Vec<String>) {
        let parts: Vec<&str> = line.split('\t').collect();
        let chrom = parts[0].to_string();
        let pos: usize = parts[1].parse().unwrap();
        let ref_base = parts[3].chars().next().unwrap();
        let alt_base = parts[4].chars().next().unwrap();
        
        // Parse INFO field for DL and DR
        let info = parts[7];
        let mut dl = 0;
        let mut dr = 0;
        for field in info.split(';') {
            if field.starts_with("DL=") {
                dl = field[3..].parse().unwrap();
            } else if field.starts_with("DR=") {
                dr = field[3..].parse().unwrap();
            }
        }
        
        // Genotypes start at column 9
        let genotypes: Vec<String> = parts[9..].iter().map(|s| s.to_string()).collect();
        
        (chrom, pos, ref_base, alt_base, dl, dr, genotypes)
    }
    
    #[test]
    fn test_vcf_biallelic_snp() {
        // 62 sites, SNP at position 31 (G/T)
        // With min_dist=30: should select the SNP (30 monomorphic on each side)
        let lines = run_vcf_test(&["biallelic_snp.fa"], 30);
        assert_eq!(lines.len(), 1, "Should find exactly 1 SNP");
        
        let (chrom, pos, ref_base, alt_base, dl, dr, genotypes) = parse_vcf_line(&lines[0]);
        assert_eq!(chrom, "biallelic_snp");
        assert_eq!(pos, 31);
        assert_eq!(ref_base, 'G');
        assert_eq!(alt_base, 'T');
        assert_eq!(dl, 30);
        assert_eq!(dr, 31);
        assert_eq!(genotypes, vec!["0", "1", "0"]);
    }
    
    #[test]
    fn test_vcf_triallelic_excluded() {
        // 62 sites, position 31 has G/T/C (triallelic)
        // Should NOT select (not biallelic)
        let lines = run_vcf_test(&["triallelic_snp.fa"], 30);
        assert_eq!(lines.len(), 0, "Triallelic site should be excluded");
    }
    
    #[test]
    fn test_vcf_two_snps_close_excluded() {
        // SNPs at positions 15 and 31 (too close for min_dist=30)
        // Neither should pass
        let lines = run_vcf_test(&["two_snps_close.fa"], 30);
        assert_eq!(lines.len(), 0, "Close SNPs should both be excluded");
    }
    
    #[test]
    fn test_vcf_two_snps_far_excluded() {
        // SNPs at positions 31 and 50 (too close for min_dist=30)
        // Neither should pass (only 18 sites between them)
        let lines = run_vcf_test(&["two_snps_far.fa"], 30);
        assert_eq!(lines.len(), 0, "SNPs with insufficient spacing should be excluded");
    }
    
    #[test]
    fn test_vcf_three_snps_small_dist() {
        // SNPs at positions 15, 31, 50 with min_dist=5
        // All three should be selected
        let lines = run_vcf_test(&["three_snps.fa"], 5);
        assert_eq!(lines.len(), 3, "All 3 SNPs should be selected with min_dist=5");
        
        // Check positions and distances
        let (_, pos1, _, _, dl1, dr1, _) = parse_vcf_line(&lines[0]);
        let (_, pos2, _, _, dl2, dr2, _) = parse_vcf_line(&lines[1]);
        let (_, pos3, _, _, dl3, dr3, _) = parse_vcf_line(&lines[2]);
        
        assert_eq!(pos1, 15);
        assert_eq!(pos2, 31);
        assert_eq!(pos3, 50);
        
        // DL/DR values
        assert_eq!(dl1, 14); // 14 sites before pos 15
        assert_eq!(dr1, 15); // 15 sites to next SNP at 31
        assert_eq!(dl2, 15); // 15 sites from SNP at 15
        assert_eq!(dr2, 18); // 18 sites to SNP at 50
        assert_eq!(dl3, 18); // 18 sites from SNP at 31
        assert_eq!(dr3, 12); // 12 sites remaining
    }
    
    #[test]
    fn test_vcf_snp_with_n_missing_genotype() {
        // SNP with N in one sample -> missing genotype but site kept
        let lines = run_vcf_test(&["snp_with_n.fa"], 30);
        assert_eq!(lines.len(), 1, "SNP with N should still be selected");
        
        let (_, _, _, _, _, _, genotypes) = parse_vcf_line(&lines[0]);
        assert_eq!(genotypes, vec!["0", "1", "."], "Sample with N should have missing genotype");
    }
    
    #[test]
    fn test_vcf_snp_with_gap_excluded() {
        // SNP with gap in one sample -> site excluded entirely
        let lines = run_vcf_test(&["snp_with_gap.fa"], 30);
        assert_eq!(lines.len(), 0, "SNP with gap should be excluded");
    }
}
