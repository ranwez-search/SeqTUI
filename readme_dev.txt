================================================================================
SeqTUI - Developer Reference Document
================================================================================

Last updated: November 2025

This document provides context for continuing development on SeqTUI, a terminal-
based FASTA alignment viewer written in Rust. It captures key design decisions,
architecture choices, and lessons learned during development.

================================================================================
PROJECT GOAL
================================================================================

SeqTUI aims to be a fast, memory-efficient terminal viewer for large FASTA 
sequence alignments, inspired by Seaview but for the terminal. Key goals:

1. Handle very large files (500MB+, millions of nucleotides per sequence)
2. Vim-style navigation for bioinformaticians comfortable with CLI
3. Support NT→AA translation with all 33 NCBI genetic codes
4. Color-coded display matching Seaview conventions
5. Minimal dependencies, easy deployment on HPC clusters

================================================================================
ARCHITECTURE
================================================================================

The codebase follows an event-driven MVC pattern:

src/
├── main.rs         - Entry point, CLI args, jemalloc allocator setup
├── lib.rs          - Module exports
├── model.rs        - Core data structures and application state
├── fasta.rs        - FASTA parsing with memory optimization
├── event.rs        - Keyboard input handling (Action enum, apply_action)
├── ui.rs           - TUI rendering with ratatui
├── controller.rs   - Main loop connecting events to state to rendering
└── genetic_code.rs - 33 NCBI genetic codes and translation logic

Key design pattern: Events -> Actions -> State mutations -> Render

================================================================================
KEY DATA STRUCTURES (model.rs)
================================================================================

Sequence {
    name: String,           // Sequence identifier
    data: Vec<u8>,          // Raw sequence bytes (NOT String - memory optimization)
}

Alignment {
    sequences: Vec<Sequence>,
    warning: Option<String>,  // Unequal length warning
}

AppState {
    alignment: Alignment,
    translated_alignment: Option<Alignment>,  // Cached AA translation
    view_mode: ViewMode,      // Nucleotide or AminoAcid
    viewport: Viewport,       // What's visible on screen
    cursor: Cursor,           // Current position
    mode: AppMode,            // Normal, Command, Search, TranslationSettings
    help_tab: HelpTab,        // Current help tab (5 tabs)
    pending_g: bool,          // For g-prefix commands
    pending_z: bool,          // For z-prefix commands (zH, zL)
    ...
}

================================================================================
PERFORMANCE OPTIMIZATIONS
================================================================================

These were critical for handling 500MB+ files (47 sequences × 11M nucleotides):

1. SEQUENCE STORAGE: Vec<u8> instead of String
   - Eliminates UTF-8 validation overhead
   - Direct byte access without bounds checking per character
   - ~30% memory reduction

2. FASTA PARSING: Bulk read for large files
   - Files >1MB: read entire file to Vec<u8>, then parse
   - Avoids per-line allocation overhead
   - Pre-allocated capacity based on file size

3. TRANSLATION: Array-based codon lookup
   - [u8; 64] array instead of HashMap for codon table
   - Inline function base_to_index() for 2-bit encoding
   - No string allocation per codon (works directly in bytes)
   - translate_sequence(&[u8]) -> Vec<u8> (no String intermediates)

4. MEMORY ALLOCATOR: jemalloc (tikv-jemallocator)
   - In theory, helps return freed memory to OS
   - In practice, the difference is minimal on modern systems
   - Kept in case it helps on some HPC/cluster environments

5. REMOVED: Parallel translation with rayon
   - Initially added for speed, but caused 15 threads overhead
   - Single-threaded is fast enough for interactive use
   - Simpler code, lower memory footprint

================================================================================
VIM NAVIGATION DESIGN
================================================================================

Navigation is designed for Vim users but also works with arrows:

ARROW-CENTRIC:
  ←↑↓→           Move one position
  Shift+←→       Half page left/right  
  Shift+↑↓       Page up/down

VIM-CENTRIC:
  h/j/k/l        Move (left/down/up/right)
  Ctrl+U/D       Half page up/down
  zH/zL          Half page left/right (Vim's horizontal scroll)
  0/$            First/last column
  g0/gm/g$       First/middle/last VISIBLE column
  <num>|         Go to column N

NOTE: Ctrl+Arrow doesn't work on macOS (captured by system for Spaces/Mission
Control), so Shift+Arrow is used instead. Both Ctrl and Shift are supported
in code for cross-platform compatibility.

================================================================================
PENDING STATE PATTERN
================================================================================

For multi-key commands (g0, gm, g$, zH, zL), we use a "pending" state:

1. User presses 'g' → set pending_g = true
2. Next key (0, m, or $) triggers the actual command
3. Command calls clear_pending() to reset

IMPORTANT: Any action triggered by a pending state must call clear_pending()
or the app becomes unresponsive (all subsequent keys go to the pending handler
which returns Action::None for unknown keys).

See: set_pending_g(), set_pending_z(), clear_pending() in model.rs

================================================================================
HELP SYSTEM
================================================================================

Tabbed help overlay with 5 sections:
- Basics: Getting started, :q, :h, :<number>
- Arrow Nav: Arrow key navigation
- Vim Nav: Vim-style navigation  
- Search: /, ?, n, N
- Translation: :asAA, :asNT, :setcode

Navigate tabs with ←/→, h/l, or Tab. Any other key closes help.

State: help_tab: HelpTab in AppState
Actions: HelpNextTab, HelpPrevTab, DismissHelp

================================================================================
TRANSLATION SYSTEM
================================================================================

Nucleotide to amino acid translation:
- 33 NCBI genetic codes supported (Standard, Vertebrate Mito, etc.)
- 3 reading frames (+1, +2, +3)
- Translation settings UI with j/k for code, h/l for frame

The translated alignment is cached in translated_alignment. When user types
:asNT, we switch view_mode back to Nucleotide and the translation stays cached.
Typing :asAA again reuses the cache unless settings changed.

Memory note: Dropping translated_alignment with jemalloc properly returns
memory to OS. Without jemalloc, memory stays allocated.

================================================================================
COLOR SCHEME
================================================================================

Nucleotides (DNA/RNA):
  A: Red background
  C: Green background
  G: Yellow background
  T/U: Blue background

Amino acids (Seaview-style, grouped by chemical property):
  Hydrophobic (AFILMVW): Yellow
  Polar (NQST): Green
  Charged+ (KRH): Magenta/Red
  Charged- (DE): Blue
  Special (CGP): Orange/Cyan/Gray
  Stop (*): Red on white

================================================================================
FILE HANDLING
================================================================================

FASTA parsing handles:
- Standard multi-line FASTA
- Lines starting with > are headers
- Sequences can span multiple lines
- Automatic uppercase conversion
- Warning if sequences have different lengths (invalid alignment)

Edge cases:
- Empty files: Error
- Files without headers: Creates "Unknown" sequence
- Very long lines: Handled (bulk read approach)

================================================================================
TESTING
================================================================================

48 unit tests covering:
- Event handling and key mappings
- FASTA parsing edge cases
- Genetic code translation
- Model state transitions
- Color assignments

Run: cargo test

================================================================================
DEPENDENCIES
================================================================================

Cargo.toml key dependencies:
- ratatui: TUI framework
- crossterm: Terminal backend (cross-platform)
- anyhow/thiserror: Error handling
- clap: CLI argument parsing
- tikv-jemallocator: Alternative memory allocator (kept for potential benefits)

================================================================================
KNOWN ISSUES / FUTURE IMPROVEMENTS
================================================================================

Potential enhancements:
1. Selection/copy to clipboard
2. Sequence statistics (GC content, length)
3. Consensus sequence display
4. Export selected region
5. Mouse support for clicking
6. Reverse complement view
7. Multiple file comparison

Performance notes:
- Initial load of 500MB file: ~2-3 seconds
- Translation of 500MB: ~1 second
- Memory usage: stable during interactive use

================================================================================
DEVELOPMENT WORKFLOW
================================================================================

# Build and run
cargo run -- test_data/alignment.fasta

# Run tests
cargo test

# Build release (for HPC)
cargo build --release

# Check for issues
cargo clippy

================================================================================
CONTACT / CONTEXT
================================================================================

This project was developed with AI assistance (Claude/Copilot). When resuming
development, provide this file as context to quickly get back up to speed on
architecture decisions and design patterns used.

Key files to review when resuming:
1. This file (readme_dev.txt)
2. src/model.rs - Core state and data structures  
3. src/event.rs - Action definitions and key handling
4. src/genetic_code.rs - Translation logic

================================================================================
