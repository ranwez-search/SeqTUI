//! # SeqTUI - Terminal Alignment Viewer
//!
//! A terminal-based viewer for FASTA sequence alignments using ratatui.
//!
//! ## Architecture
//!
//! The application follows an event-driven architecture with clear separation:
//! - `model`: Data structures for sequences, viewport, and application state
//! - `fasta`: FASTA file parsing and validation
//! - `event`: Keyboard event handling (Vim-style navigation)
//! - `ui`: TUI rendering with ratatui
//! - `controller`: Orchestration of state transitions
//!
//! ## Future Extensions
//!
//! The architecture is designed to support:
//! - Amino acid color schemes
//! - Large alignment handling
//! - Pattern search
//! - Column navigation via `:number`
//! - File browser panel
//! - Sequence/site filtering
//! - Codon coloring and NT to AA translation
//! - Export current view as FASTA

pub mod controller;
pub mod event;
pub mod fasta;
pub mod model;
pub mod ui;
