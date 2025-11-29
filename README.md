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
git clone https://github.com/your-username/SeqTUI.git
cd SeqTUI

# Build and install
cargo install --path .
```

## Usage

```bash
# View an alignment (format auto-detected from extension)
seqtui sequences.fasta
seqtui alignment.phy
seqtui data.nex

# Force a specific format
seqtui -f nexus myfile.txt
seqtui -f phylip alignment.dat

# Or run with cargo
cargo run -- sequences.fasta
```

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
