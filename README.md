# SeqTUI - Terminal Alignment Viewer

A terminal-based viewer for FASTA sequence alignments written in Rust using [ratatui](https://ratatui.rs/).

## Features

- 🧬 **FASTA Support**: Load and visualize aligned DNA/protein sequences
- 🎨 **Color Coded**: Nucleotides displayed with distinct background colors
  - A: Red
  - C: Green
  - G: Yellow
  - T: Blue
  - Others: Gray
- 📜 **Sticky Names**: Sequence identifiers remain visible while scrolling horizontally
- ⌨️ **Vim-style Navigation**: Intuitive keyboard controls
- 🔍 **Auto-centering**: Column stays centered when scrolling

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

| Key | Action |
|-----|--------|
| `i` | Move up |
| `k` | Move down |
| `j` | Move right (next column) |
| `l` | Move left (previous column) |
| `:q` | Quit |
| `:123` | Jump to column 123 |
| `Ctrl+C` | Emergency quit |

Arrow keys also work for navigation.

## Architecture

The application follows an event-driven architecture:

```
src/
├── main.rs        # Entry point and CLI argument parsing
├── lib.rs         # Module exports
├── model.rs       # Data structures (Sequence, Alignment, Viewport, AppState)
├── fasta.rs       # FASTA file parsing
├── event.rs       # Keyboard event handling
├── ui.rs          # TUI rendering with ratatui
└── controller.rs  # Main application loop
```

## Future Extensions

The architecture is designed to support future enhancements:

- [ ] Amino acid color schemes
- [ ] Large alignment optimization
- [ ] Pattern search (/)
- [ ] File browser panel
- [ ] Sequence/site filtering
- [ ] Codon coloring
- [ ] NT → AA translation
- [ ] Export current view as FASTA

## Development

```bash
# Run tests
cargo test

# Run with test data
cargo run -- test_data/alignment.fasta

# Build release version
cargo build --release
```

## License

MIT
