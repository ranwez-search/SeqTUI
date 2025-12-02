//! # SeqTUI - Terminal Alignment Viewer
//!
//! A terminal-based viewer for FASTA sequence alignments using ratatui.
//!
//! ## Architecture
//!
//! The application follows an event-driven architecture with clear separation:
//! - `model`: Data structures for sequences, viewport, and application state
//! - `fasta`: FASTA file parsing and validation
//! - `formats`: Multi-format parsing (FASTA, PHYLIP, NEXUS)
//! - `event`: Keyboard event handling (Vim-style navigation)
//! - `ui`: TUI rendering with ratatui
//! - `controller`: Orchestration of state transitions
//! - `genetic_code`: NCBI genetic codes and translation logic
//!
//! ## Supported File Formats
//!
//! - FASTA (.fasta, .fa, .fna, .faa)
//! - PHYLIP (.phy, .phylip) - sequential and interleaved
//! - NEXUS (.nex, .nexus, .nxs)

pub mod controller;
pub mod event;
pub mod formats;
pub mod genetic_code;
pub mod model;
pub mod ui;
