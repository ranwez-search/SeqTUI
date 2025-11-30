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
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Sequence file to display (FASTA, PHYLIP, or NEXUS format)
    file: PathBuf,

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

    /// Reading frame for translation (1-6, default: 1)
    /// Frames 1-3 are forward (+1, +2, +3), 4-6 are reverse complement
    #[arg(short = 'r', long = "reading-frame", default_value = "1")]
    reading_frame: u8,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let file_path = args.file;
    let forced_format: Option<FileFormat> = args.format.into();

    // Validate reading frame (1-6)
    if args.reading_frame < 1 || args.reading_frame > 6 {
        anyhow::bail!("Reading frame must be 1-6 (got {})", args.reading_frame);
    }

    // Validate genetic code (1-33, with some gaps)
    if args.genetic_code < 1 || args.genetic_code > 33 {
        anyhow::bail!("Genetic code must be 1-33 (got {})", args.genetic_code);
    }

    // CLI mode: output to file/stdout
    if let Some(output) = args.output {
        run_cli_mode(&file_path, forced_format, &output, args.translate, args.genetic_code, args.reading_frame)?;
    } else {
        // TUI mode with optional preset translation settings
        run_app_with_loading(
            file_path,
            forced_format,
            if args.translate || args.genetic_code != 1 || args.reading_frame != 1 {
                Some((args.genetic_code, args.reading_frame))
            } else {
                None
            },
        )?;
    }

    Ok(())
}
