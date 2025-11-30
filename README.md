# SeqTUI - Terminal Alignment Viewer

A terminal-based viewer for sequence alignments written in Rust using [ratatui](https://ratatui.rs/).

## Features

- 🧬 **Multi-Format Support**: FASTA, PHYLIP, and NEXUS formats with auto-detection
- 🎨 **Color Coded**: Nucleotides and amino acids displayed with distinct background colors
  - **Nucleotides**: A (Red), C (Green), G (Yellow), T/U (Blue)
  - **Amino Acids**: Seaview-style coloring by chemical properties
- 🔄 **NT → AA Translation**: Translate nucleotide sequences to amino acids with 33 genetic codes
- 📜 **Sticky Names**: Sequence identifiers remain visible while scrolling horizontally
- ⌨️ **Vim-style Navigation**: Full Vim-like controls (h/j/k/l, w/b/e, Ctrl+U/D, g0/gm/g$)
- 🖱️ **Arrow Navigation**: Shift+arrows for page/half-page scrolling
- 🔍 **Pattern Search**: Search forward (`/`) and backward (`?`) in sequences and names
- ❓ **Tabbed Help**: Press `:h` for organized help with 5 sections
- 🚀 **Large File Support**: Tested on 500MB+ alignments

## Installation

```bash
# Clone the repository
git clone https://github.com/ranwez-search/SeqTUI.git
cd SeqTUI

# Build and install
cargo install --path .
```

## Usage

SeqTUI works in two modes:
- **TUI mode** (default): Interactive terminal viewer
- **CLI mode** (with `-o`): Command-line converter for batch processing

### TUI Mode (Interactive Viewer)

```bash
# View an alignment (format auto-detected from extension)
seqtui sequences.fasta
seqtui alignment.phy
seqtui data.nex

# Force a specific format
seqtui -f nexus myfile.txt
seqtui -f phylip alignment.dat

# Preset translation settings for TUI
seqtui sequences.fasta -g 2 -r 1    # Open with genetic code 2 preset
```

### CLI Mode (Batch Processing)

Use `-o` to output to a file (or `-` for stdout) instead of opening the TUI.
**Single-line FASTA output** makes sequences easy to process with standard Unix tools:

```bash
# Convert to single-line FASTA
seqtui sequences.fasta -o sequences_1L.fasta

# Translate to amino acids
seqtui sequences.fasta -t -o sequences_AA.fasta

# Translate with specific genetic code and reading frame
seqtui sequences.fasta -t -g 2 -r 1 -o sequences_AA.fasta
```

**Unix pipeline examples** (single-line FASTA makes this easy!):

```bash
# Check for internal stop codons in coding sequences
seqtui sequences.fasta -t -o - | grep "\*."

# Create a small test set (first 10 sequences, 500 bp each)
seqtui large_alignment.nex -o - | head -20 | cut -c1-500 > test_seq.fasta

# Extract subset of sequences by ID
seqtui sequences.nex -o - | grep -A1 -f seq_ids.txt > subset.fasta

# Count sequences
seqtui alignment.phy -o - | grep -c "^>"

# Batch translate all files in a directory
for f in *.fasta; do seqtui "$f" -t -o "${f%.fasta}_AA.fasta"; done
```

### CLI Options

| Option | Long | Description |
|--------|------|-------------|
| `-o` | `--output` | Output file (use `-` for stdout). Enables CLI mode |
| `-t` | `--translate` | Translate nucleotides to amino acids |
| `-g` | `--genetic-code` | Genetic code (1-33, default: 1 = Standard) |
| `-r` | `--reading-frame` | Reading frame (1-3, default: 1) |
| `-f` | `--format` | Force input format (fasta, phylip, nexus, auto) |

## Supported Formats

| Format | Extensions | Features |
|--------|------------|----------|
| **FASTA** | `.fasta`, `.fa`, `.fna`, `.faa`, `.fas` | Multi-line sequences |
| **PHYLIP** | `.phy`, `.phylip` | Sequential and interleaved |
| **NEXUS** | `.nex`, `.nexus`, `.nxs` | DATA/CHARACTERS blocks, MATCHCHAR support |

## Navigation

### Arrow Keys
| Key | Action |
|-----|--------|
| `←↑↓→` | Move one position |
| `Shift+←→` | Half page left/right |
| `Shift+↑↓` | Full page up/down |
| `Home` / `End` | First / last column |
| `PgUp` / `PgDn` | Full page up/down |

### Vim-style
| Key | Action |
|-----|--------|
| `h` / `j` / `k` / `l` | Move left/down/up/right |
| `Ctrl+U` / `Ctrl+D` | Half page up/down |
| `zH` / `zL` | Half page left/right |
| `0` / `$` | First / last column |
| `g0` / `gm` / `g$` | First/middle/last visible column |
| `<num>\|` | Go to column (e.g., `50\|`) |
| `w` / `b` / `e` | Next/previous/end of word |

## Search

| Key | Action |
|-----|--------|
| `/pattern` | Search forward |
| `?pattern` | Search backward |
| `n` | Next match |
| `N` | Previous match |

## Commands

| Command | Action |
|---------|--------|
| `:q` | Quit |
| `:h` | Toggle help overlay |
| `:<number>` | Go to sequence/row |
| `:w file.fa` | Save current view to FASTA (single-line format) |
| `:asAA` | Translate nucleotides to amino acids |
| `:asNT` | Switch back to nucleotide view |
| `:setcode` | Change genetic code and reading frame |

## Saving Alignments

### In TUI Mode
Use `:w filename.fasta` to save the current view:
- Saves NT or AA sequences depending on current view mode
- Sequences are written on single lines (convenient for bash processing)
- Example: `:w Loc256_AA.fasta`

### In CLI Mode
Use `-o` for batch conversion — see [CLI Mode](#cli-mode-batch-processing) for examples.

## Translation

For nucleotide alignments, use `:asAA` to translate to amino acids:
- Choose from 33 NCBI genetic codes
- Select reading frame (+1, +2, +3)
- Use `j`/`k` to browse codes, `h`/`l` to change frame
- Press `Enter` to translate, `Esc` to cancel

Use `:asNT` to switch back to the nucleotide view.

## Architecture

The application follows an event-driven MVC architecture:

```
src/
├── main.rs         # Entry point and CLI argument parsing
├── lib.rs          # Module exports
├── model.rs        # Data structures (Sequence, Alignment, Viewport, AppState)
├── fasta.rs        # FASTA file parsing
├── formats/        # Multi-format support
│   ├── mod.rs      # Format detection and unified parsing
│   ├── nexus.rs    # NEXUS format parser (token-based)
│   └── phylip.rs   # PHYLIP format parser
├── event.rs        # Keyboard event handling
├── ui.rs           # TUI rendering with ratatui
├── controller.rs   # Main application loop
└── genetic_code.rs # Genetic code tables and translation
```

## Development

```bash
# Run tests (74+ tests covering all formats)
cargo test

# Run with test data
cargo run -- test_data/alignment.fasta
cargo run -- test_data/LOC_01790.nex

# Build release version
cargo build --release
```

## Building for Linux (HPC)

```bash
# On Linux or cross-compile
cargo build --release --target x86_64-unknown-linux-gnu
```

## Note

⚠️ This is a side project developed with extensive use of AI assistants (GitHub Copilot and Claude Opus 4.5). It was my first "vibe coding" experiment. While functional, please use with appropriate caution and feel free to report any issues.

## Related Projects & Inspiration

SeqTUI was inspired by several great tools:

- **[Seaview](https://doua.prabi.fr/software/seaview)** - The reference GUI alignment viewer. I'm a huge fan and it has many more features. SeqTUI was born from needing a quick alignment view when working on HPC clusters via terminal.
- **[Vim](https://www.vim.org/)** - The navigation philosophy and keybindings
- **[ratatui](https://ratatui.rs/)** - The excellent Rust TUI framework powering this app
- **[Awesome ratatui](https://github.com/ratatui/awesome-ratatui)** - Inspiration from the ecosystem, especially:
  - [termscp](https://github.com/veeso/termscp) - Terminal file transfer
  - [tgv](https://github.com/chroma-pipeline/tgv) - Terminal genome viewer

## License

MIT
