# SeqTUI - Terminal Alignment Viewer

A terminal-based viewer for FASTA sequence alignments written in Rust using [ratatui](https://ratatui.rs/).

## Features

- 🧬 **FASTA Support**: Load and visualize aligned DNA/RNA and protein sequences
- 🎨 **Color Coded**: Nucleotides and amino acids displayed with distinct background colors
  - **Nucleotides**: A (Red), C (Green), G (Yellow), T/U (Blue)
  - **Amino Acids**: Seaview-style coloring by chemical properties
- 🔄 **NT → AA Translation**: Translate nucleotide sequences to amino acids with 33 genetic codes
- 📜 **Sticky Names**: Sequence identifiers remain visible while scrolling horizontally
- ⌨️ **Vim-style Navigation**: Full Vim-like controls (h/j/k/l, Ctrl+U/D, zH/zL, g0/gm/g$)
- 🖱️ **Arrow Navigation**: Shift+arrows for page/half-page scrolling
- 🔍 **Pattern Search**: Search forward (`/`) and backward (`?`) in sequences and names
- ❓ **Tabbed Help**: Press `:h` for organized help with 5 sections
- 🚀 **Large File Support**: Optimized for 500MB+ alignments with jemalloc allocator

## Installation

```bash
# Clone the repository
git clone https://github.com/your-username/SeqTUI.git
cd SeqTUI

# Build and install
cargo install --path .
```

## Usage

```bash
# View an alignment
seqtui sequences.fasta

# Or run with cargo
cargo run -- sequences.fasta
```

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
| `:asAA` | Translate nucleotides to amino acids |
| `:asNT` | Switch back to nucleotide view |
| `:setcode` | Change genetic code and reading frame |

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
├── event.rs        # Keyboard event handling
├── ui.rs           # TUI rendering with ratatui
├── controller.rs   # Main application loop
└── genetic_code.rs # Genetic code tables and translation
```

## Development

```bash
# Run tests
cargo test

# Run with test data
cargo run -- test_data/alignment.fasta

# Build release version
cargo build --release
```

## Building for Linux (HPC)

```bash
# On Linux or cross-compile
cargo build --release --target x86_64-unknown-linux-gnu
```

## License

MIT
