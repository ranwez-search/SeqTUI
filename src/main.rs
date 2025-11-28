//! SeqTUI - Terminal Alignment Viewer
//!
//! A terminal-based viewer for FASTA sequence alignments.
//!
//! ## Usage
//!
//! ```bash
//! seqtui <fasta_file>
//! ```
//!
//! ## Navigation (Vim-style)
//!
//! - `i`: Move up
//! - `k`: Move down
//! - `j`: Move right (next column)
//! - `l`: Move left (previous column)
//! - `:q`: Quit
//! - `:123`: Jump to column 123

use anyhow::{Context, Result};
use clap::Parser;
use std::path::PathBuf;

use seqtui::controller::run_app;
use seqtui::fasta::parse_fasta_file;
use seqtui::model::AppState;

/// SeqTUI - Terminal Alignment Viewer
///
/// A terminal-based viewer for FASTA sequence alignments with Vim-style navigation
/// and nucleotide color coding.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// FASTA file(s) to display
    #[arg(required = true)]
    files: Vec<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // For now, we only support single file
    // Future: support multiple files or directory browsing
    let file_path = &args.files[0];

    // Parse the FASTA file
    let alignment = parse_fasta_file(file_path)
        .with_context(|| format!("Failed to parse FASTA file: {}", file_path.display()))?;

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
