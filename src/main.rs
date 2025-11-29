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

use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use std::path::PathBuf;

use seqtui::controller::run_app;
use seqtui::formats::{parse_file_with_options, FileFormat};
use seqtui::model::AppState;

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
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Sequence file to display (FASTA, PHYLIP, or NEXUS format)
    file: PathBuf,

    /// Force a specific file format (overrides auto-detection)
    #[arg(short = 'f', long = "format", value_enum, default_value = "auto")]
    format: FormatArg,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let file_path = &args.file;
    let forced_format: Option<FileFormat> = args.format.into();

    // Parse the sequence file
    let alignment = parse_file_with_options(file_path, forced_format)
        .with_context(|| format!("Failed to parse file: {}", file_path.display()))?;

    // Print info before entering TUI (will be hidden)
    if !alignment.is_valid_alignment {
        eprintln!(
            "Warning: {}",
            alignment.warning.as_ref().unwrap_or(&String::new())
        );
    }

    // Extract file name (basename without extension)
    let file_name = file_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("alignment")
        .to_string();

    // Create application state
    let state = AppState::new(alignment, file_name);

    // Run the TUI application
    run_app(state)?;

    Ok(())
}
