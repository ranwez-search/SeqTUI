================================================================================
SeqTUI - Developer Reference Document
================================================================================

Last updated: November 2025

This document provides context for continuing development on SeqTUI, a terminal-
based sequence alignment viewer (FASTA, PHYLIP, NEXUS) written in Rust. It 
captures key design decisions, architecture choices, and lessons learned.

================================================================================
PROJECT GOAL
================================================================================

SeqTUI aims to be a fast, memory-efficient terminal viewer for large 
sequence alignments (FASTA, PHYLIP, NEXUS), inspired by Seaview but for 
the terminal. Key goals:

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
├── formats/        - Multi-format support module
│   ├── mod.rs      - Format detection (extension + content) and unified API
│   ├── nexus.rs    - NEXUS parser (token-based per spec)
│   └── phylip.rs   - PHYLIP parser (sequential + interleaved)
├── event.rs        - Keyboard input handling (Action enum, apply_action)
├── ui.rs           - TUI rendering with ratatui
├── controller.rs   - Main loop, background loading, channel-based messaging
└── genetic_code.rs - 33 NCBI genetic codes and translation logic

Key design pattern: Events -> Actions -> State mutations -> Render
Async pattern: Background thread -> Channel -> Main loop polls -> State update

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
    cached_translation_code_id: Option<u8>,   // Code used for cached translation
    cached_translation_frame: Option<usize>,  // Frame used for cached translation
    view_mode: ViewMode,      // Nucleotide or AminoAcid
    loading_state: LoadingState,  // Ready, LoadingFile, Translating
    spinner_frame: usize,     // Animation frame (0-3)
    viewport: Viewport,       // What's visible on screen
    cursor: Cursor,           // Current position
    mode: AppMode,            // Normal, Command, Search, TranslationSettings
    help_tab: HelpTab,        // Current help tab (5 tabs)
    pending_g: bool,          // For g-prefix commands
    pending_z: bool,          // For z-prefix commands (zH, zL)
    ...
}

LoadingState {
    Ready,                    // No loading in progress
    LoadingFile { path, message, sequences_loaded },
    Translating { message, sequences_done, total },
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
ASYNC LOADING ARCHITECTURE
================================================================================

The TUI opens IMMEDIATELY when the user runs seqtui. File parsing happens in
a background thread while a loading spinner is displayed.

Components:

1. LoadingState enum (model.rs)
   - Ready: No loading, alignment is displayed
   - LoadingFile: Parsing in progress, shows spinner
   - Translating: Translation in progress (future use)

2. LoadMessage enum (controller.rs)
   - Complete(Alignment): Parsing succeeded
   - Error(String): Parsing failed
   - Progress { sequences_loaded }: Streaming updates (future use)

3. Background loading flow:
   a. main.rs calls run_app_with_loading(path, format)
   b. Controller creates AppState in LoadingFile state
   c. Controller spawns std::thread for parsing
   d. TUI renders immediately with spinner overlay
   e. Main loop polls channel (non-blocking try_recv)
   f. On LoadMessage::Complete, state.set_alignment() is called
   g. Spinner disappears, alignment is shown

4. Spinner animation:
   - Braille characters: ⠋ ⠙ ⠹ ⠸ (4 frames)
   - tick_spinner() advances frame each render loop (~50ms)
   - spinner_char() returns current frame character

5. Error handling:
   - If parsing fails, LoadMessage::Error is sent
   - state.set_loading_error() displays error in status bar
   - User can quit with :q

Future: Could add streaming parser that sends Progress messages as sequences
are parsed, updating the count in the loading overlay.

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

TRANSLATION CACHING:
The translated alignment is cached in `translated_alignment` along with 
metadata tracking which settings were used:
- cached_translation_code_id: Option<u8> - Genetic code ID used
- cached_translation_frame: Option<usize> - Frame used (0, 1, or 2)

When user types :asNT, we switch view_mode back to Nucleotide but KEEP the
cached translation. When typing :asAA again:
1. has_valid_cached_translation() checks if cached settings match current
2. If match: switch_to_cached_aa_view() - instant, no recomputation
3. If no match: start_background_translation() - recompute in background

This enables rapid NT↔AA toggling without recomputation, which is important
for large alignments where translation can take several seconds.

CACHE INVALIDATION:
Cache is invalidated when genetic_code_id or frame changes:
- User opens :setcode dialog and changes settings
- Cache metadata won't match, triggering recomputation

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

Multi-format support with auto-detection:

FORMAT DETECTION STRATEGY (in parse_file_with_options):
-------------------------------------------------------

The detection follows a cascading fallback strategy:

1. EXPLICIT FORMAT (-f/--format option)
   If user specifies -f nexus, we use NEXUS parser directly.
   If parsing fails, we return the error (no fallback).

2. FILE EXTENSION (if no -f option)
   We try the parser matching the extension (.fasta → FASTA parser).
   IMPORTANT: If extension-based parsing FAILS, we SILENTLY fall through
   to content detection. No warning is displayed.
   
   Example: seq.fasta containing NEXUS data
   - Try FASTA parser (fails because #NEXUS is not a valid FASTA header)
   - Fall through to step 3

3. CONTENT-BASED DETECTION
   Examine first non-empty line:
   - Starts with "#NEXUS" (case-insensitive) → NEXUS
   - Starts with ">" → FASTA  
   - Two integers (ntax nchar) → PHYLIP
   
   If detected, parse with that format.
   If parsing fails, return the error (no further fallback).

4. TRY ALL PARSERS (last resort)
   If content detection returns None, try each parser in order:
   FASTA → NEXUS → PHYLIP
   Return first success, or UnknownFormat error if all fail.

CURRENT BEHAVIOR SUMMARY:
- Extension is a HINT, not authoritative
- Extension mismatch produces NO warning (silent fallback)
- Content signature is trusted when found
- User can always override with -f/--format

POTENTIAL IMPROVEMENT:
Could add warning when extension doesn't match detected/successful format:
  "Warning: file.fasta was parsed as NEXUS (extension suggests FASTA)"
Currently NOT implemented - parsing succeeds silently.

Format detection priority:
1. Explicit -f/--format CLI option (fasta, phylip, nexus, auto)
2. File extension (.fasta, .phy, .nex, etc.) - with silent fallback on failure
3. Content detection (looks for format signatures)
4. Try all parsers as fallback

FASTA parsing handles:
- Standard multi-line FASTA
- Lines starting with > are headers
- Sequences can span multiple lines
- Automatic uppercase conversion
- Warning if sequences have different lengths (invalid alignment)

PHYLIP parsing handles:
- Sequential format (all of sequence on consecutive lines)
- Interleaved format (detected by line count vs NCHAR)
- Relaxed names (any length, not just 10 chars)
- Strict 10-char names for legacy files

NEXUS parsing handles:
- Token-based parsing per NEXUS specification
- DATA and CHARACTERS blocks
- Sequential and INTERLEAVE formats
- MATCHCHAR substitution (e.g., '.' = same as reference sequence)
- Multi-line FORMAT commands
- Inline comments like [1], [annotation]
- Case-insensitive commands
- Quoted sequence names

Edge cases:
- Empty files: Error
- Files without headers: Creates "Unknown" sequence (FASTA)
- Very long lines: Handled (bulk read approach)
- Unknown format: Helpful error with format hints

================================================================================
TESTING
================================================================================

75+ unit tests covering:
- Event handling and key mappings
- FASTA parsing edge cases
- PHYLIP parsing (sequential, interleaved, relaxed names)
- NEXUS parsing (simple, interleaved, quoted names, MATCHCHAR)
- Format detection (extension, content, fallback)
- Real-world file tests (LOC_01790.nex with 27 sequences)
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
2. src/model.rs - Core state and data structures (incl. LoadingState)
3. src/controller.rs - Main loop, background loading, LoadMessage channel
4. src/event.rs - Action definitions and key handling
5. src/genetic_code.rs - Translation logic
6. src/formats/mod.rs - Format detection and unified parsing API
7. src/formats/nexus.rs - Token-based NEXUS parser

================================================================================
