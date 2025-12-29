# SeqTUI - Terminal Alignment Viewer & Toolkit

A fast terminal-based viewer and command-line toolkit for sequences. View, translate, convert (to FASTA), and combine sequences aligned or not — all from the terminal.

Written in Rust using [ratatui](https://ratatui.rs/).

**Author:** [Vincent Ranwez](https://github.com/ranwez)

<p align="center">
  <img src="docs/screenshot-help-translation.jpg"
       alt="SeqTUI help and translation view"
       width="80%">
</p>

## Features

- 🧬 **Multi-Format Support**: FASTA, PHYLIP, and NEXUS formats with auto-detection
- 🎨 **Color-Coded**: Nucleotides and amino acids displayed with distinct background colors
- 🔄 **NT → AA Translation**: 33 NCBI genetic codes, handling ambiguity codes (R, Y, N)
- 🧩 **Concatenation & Supermatrix**: Combine multiple alignments, fill missing with gaps
- 🧬 **SNP Extraction**: Export isolated biallelic SNPs to VCF with flanking distance filter
- 📜 **Sticky Names**: Sequence identifiers remain visible while scrolling
- ⌨️ **Vim-style Navigation**: h/j/k/l, w/b/e, Ctrl+U/D, search with `/` and `?`
- 🔧 **Command-Line Mode**: Batch convert, translate, combine — pipe-friendly output
- 🚀 **Large File Support**: Tested on alignments larger than 500 MB 

## Quick start

```bash
# View an alignment interactively
seqtui LOC_01790.nex

# Translate to amino acids
seqtui alignment.fasta -t -o out_AA.fasta
```

## Download binary releases or build from source

Precompiled binaries are available from the **Releases** section of the GitHub repository:
https://github.com/ranwez-search/SeqTUI/releases

**macOS note:** on first launch, macOS may block the application. In this case, allow it via  
*System Settings → Privacy & Security*, then run it again.

**terminal requirements:** SeqTUI requires a terminal that supports ANSI colors and cursor movement (most modern terminals do). By default, SeqTUI uses an ASCII-safe interface for maximum compatibility. Optional Unicode glyphs can be enabled using the `--fancy` option. While Unicode glyphs provide a nicer visual rendering, they may cause display issues on some terminal configurations, especially on Windows.

SeqTUI has been tested on macOS using [WezTerm](https://wezterm.org), including remote use over SSH to Linux systems, and on Windows 11 using [Windows Terminal](https://learn.microsoft.com/windows/terminal/).

Alternatively, SeqTUI can be built from source. First install Rust and Cargo by following the official instructions at  
https://www.rust-lang.org/tools/install

Then clone the repository and build SeqTUI locally:

```bash
git clone https://github.com/ranwez-search/SeqTUI.git
cd SeqTUI
cargo install --path .
```

## Usage

SeqTUI works in two modes:
- **Interactive Viewer** (default): Opens a full-screen terminal interface
- **Command-Line Mode** (with `-o`): Batch processing for scripts and pipelines

Small test data files are available in the `test_data/` directory of the repository.

---

## Command-Line Mode

Use `-o` to output to a file (or `-` for stdout) instead of opening the interactive viewer.
**Single-line FASTA output** makes sequences easy to process with standard Unix tools.

### Basic Examples

```bash
# Convert to single-line FASTA
seqtui sequences.fasta -o sequences_1L.fasta

# Translate to amino acids
seqtui sequences.fasta -t -o sequences_AA.fasta

# Translate with specific genetic code and reading frame
seqtui sequences.fasta -t -g 2 -r 1 -o sequences_AA.fasta
```

### Unix Pipeline Examples

```bash
# Check for internal stop codons in coding sequences
seqtui sequences.fasta -t -o - | grep "\*."

# Create a small test set (first 10 sequences, 500 bp each)
seqtui large_alignment.nex -o - | head -20 | cut -c1-500 > test_seq.fasta

# Extract subset of sequences by ID
seqtui sequences.nex -o - | grep -A1 -w -f seq_ids.txt > subset.fasta

# Count sequences
seqtui alignment.phy -o - | grep -c "^>"

# Batch translate all files in a directory
for f in *.fasta; do seqtui "$f" -t -o "${f%.fasta}_AA.fasta"; done
```

### Sequence Concatenation & Supermatrix

Concatenate multiple alignment files by matching sequence IDs:

```bash
# Basic concatenation (sequences matched by ID)
seqtui gene1.fasta gene2.fasta gene3.fasta -o concatenated.fasta

# Mix different formats (FASTA, PHYLIP, NEXUS) - auto-detected!
seqtui gene1.fasta gene2.phy gene3.nex -o concatenated.fasta

# Supermatrix: fill missing sequences with gaps
seqtui gene*.fasta -s -o supermatrix.fasta

# With partition file for phylogenetic analysis
seqtui gene*.fasta -s -p partitions.txt -o supermatrix.fasta
# Creates: partitions.txt with "gene1 = 1-500", "gene2 = 501-1200", etc.

# Translate all genes and build AA supermatrix
seqtui gene*.fasta -s -t -o supermatrix_AA.fasta
```

**Sequence ID matching with delimiter** — when sequence names have prefixes/suffixes:

```bash
# Files have: Human_ENS001, Human_LOC789, Mouse_ENS002, Mouse_LOC456...
# Match on species name (first field before "_")
seqtui gene1.fasta gene2.fasta -d "_" -s -o supermatrix.fasta
# Output sequences: Human, Mouse, ... (matched across files)

# Keep multiple fields: Ae_bicornis_contig257 → Ae_bicornis
seqtui gene*.fasta -f 1,2 -d "_" -s -o supermatrix.fasta
```

### SNP Extraction (VCF Export)

Extract isolated biallelic SNPs from alignments using a minimum flanking monomorphic distance filter:

```bash
# Extract SNPs with at least 300 monomorphic sites on each side
seqtui alignment.fasta -v 300 -o snps.vcf

# Process multiple alignments (each becomes a separate CHROM)
seqtui gene*.fasta -v 100 -o all_snps.vcf

# With sequence ID matching via delimiter
seqtui gene*.fasta -v 100 -d "_" -o snps.vcf
```

**VCF output features:**
- Reference = first sequence of first file
- Samples sorted alphabetically (reference first)
- Haploid genotypes: `0` (ref), `1` (alt), `.` (missing)
- `DL` and `DR` in INFO: distance to nearest polymorphic site (for filtering)
- Only isolated biallelic SNPs are exported (polymorphic sites reduce DL/DR)
- Sites with gaps are excluded; N/? become missing genotypes

### CLI Options

| Option | Long | Description |
|--------|------|-------------|
| `-o` | `--output` | Output file in sorted FASTA (triggers CLI mode). Use `-` for stdout |
| | `--format` | Force input format (fasta, phylip, nexus). Default: auto-detect |
| | `--force` | Proceed despite warnings (ID mismatches, suspect sequences) |
| `-d` | `--delimiter` | Delimiter for splitting sequence IDs (default: `_` when -f is used) |
| `-f` | `--fields` | Fields to keep from IDs (1-based, comma-separated). Ex: `-f 1,2` |
| `-s` | `--supermatrix` | Fill missing sequences with a character (default: `-`) |
| `-p` | `--partitions` | Write partition file for phylogenetic analysis |
| `-t` | `--translate` | Translate nucleotides to amino acids |
| `-g` | `--genetic-code` | Genetic code (1-33, default: 1 = Standard) |
| `-r` | `--reading-frame` | Reading frame (1-3, default: 1) |
| `-v` | `--vcf` | Extract isolated biallelic SNPs to VCF (value = min flanking distance) |
| | `--fancy` | Enable fancy Unicode glyphs in the interactive TUI (may cause issues on some terminals, especially on Windows) |

### Fancy UI shortcut

On Unix-like systems (Linux, macOS), you may prefer to always use the fancy Unicode interface. In this case, you can define a shell alias:

```bash
alias seqtui='seqtui --fancy'
```

This keeps the default behavior safe and portable while allowing a richer interface when explicitly requested.

---

## Interactive Viewer

By default, SeqTUI opens a full-screen terminal interface for exploring alignments.

### Opening Files

```bash
# Launch file browser (no arguments)
seqtui

# View an alignment (format auto-detected from extension)
seqtui sequences.fasta
seqtui alignment.phy
seqtui data.nex

# Force a specific format
seqtui --format nexus myfile.txt
seqtui --format phylip alignment.dat

# Preset translation settings
seqtui sequences.fasta -g 2 -r 1    # Open with genetic code 2 preset
```

### Navigation

#### Arrow Keys
| Key | Action |
|-----|--------|
| `←↑↓→` | Move one position |
| `Shift+←→` | Half page left/right |
| `Shift+↑↓` | Full page up/down |
| `Home` / `End` | First / last column |
| `PgUp` / `PgDn` | Full page up/down |

#### Vim-style
| Key | Action |
|-----|--------|
| `h` / `j` / `k` / `l` | Move left/down/up/right |
| `Ctrl+U` / `Ctrl+D` | Half page up/down |
| `zH` / `zL` | Half page left/right |
| `0` / `$` | First / last column |
| `g0` / `gm` / `g$` | First/middle/last visible column |
| `<num>\|` | Go to column (e.g., `50\|`) |
| `w` / `b` / `e` | Next/previous/end of word |

### Search

| Key | Action |
|-----|--------|
| `/pattern` | Search forward |
| `?pattern` | Search backward |
| `n` | Next match |
| `N` | Previous match |

### Commands

| Command | Action |
|---------|--------|
| `:q` | Quit |
| `:h` | Toggle help overlay |
| `:<number>` | Go to sequence/row |
| `:e` | Open file browser |
| `:w file.fa` | Save current view to FASTA |
| `:asAA` | Translate nucleotides to amino acids |
| `:asNT` | Switch back to nucleotide view |
| `:setcode` | Change genetic code and reading frame |

### Translation

Use `:asAA` to translate nucleotides to amino acids:
- Choose from 33 NCBI genetic codes
- Select reading frame (+1, +2, +3)
- Use `j`/`k` to browse codes, `h`/`l` to change frame
- Press `Enter` to translate, `Esc` to cancel

Use `:asNT` to switch back to the nucleotide view.

### Saving

Use `:w filename.fasta` to save the current view:
- Saves NT or AA sequences depending on current view mode
- Sequences are written on single lines (convenient for bash processing)
- Example: `:w Loc256_AA.fasta`

---

## Supported Formats

| Format | Extensions | Features |
|--------|------------|----------|
| **FASTA** | `.fasta`, `.fa`, `.fna`, `.faa`, `.fas` | Multi-line sequences |
| **PHYLIP** | `.phy`, `.phylip` | Sequential and interleaved |
| **NEXUS** | `.nex`, `.nexus`, `.nxs` | DATA/CHARACTERS blocks, MATCHCHAR support |

## Architecture

The application follows an event-driven MVC architecture. The internal architecture is shown below for developers interested in extending or modifying SeqTUI:

```
src/
├── main.rs         # Entry point and CLI argument parsing
├── lib.rs          # Module exports
├── model.rs        # Data structures (Sequence, Alignment, Viewport, AppState)
├── formats/        # Multi-format support
│   ├── mod.rs      # Format detection and unified parsing
│   ├── fasta.rs    # FASTA format parser
│   ├── nexus.rs    # NEXUS format parser (token-based)
│   └── phylip.rs   # PHYLIP format parser
├── event.rs        # Keyboard event handling
├── ui.rs           # TUI rendering with ratatui
├── controller.rs   # Main application loop
└── genetic_code.rs # Genetic code tables and translation
```

## Development

```bash
# Run tests (96 tests covering all formats and SNP extraction functionality)
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

⚠️ This is a side project developed with extensive use of AI assistants (Claude Opus 4.5 via GitHub Copilot) as a first “vibe coding” experiment.
SeqTUI has been thoroughly tested on real data with a strong focus on functionality and user experience. Please report any issues or suggestions.

## Related Projects & Inspiration

SeqTUI was inspired by several great tools:

- **[Seaview](https://doua.prabi.fr/software/seaview)** - The reference GUI alignment viewer. I'm a huge fan and it has many more features. SeqTUI was born from needing a quick alignment view when working on HPC clusters via terminal.
- **[Vim](https://www.vim.org/)** - The navigation philosophy and keybindings
- **[ratatui](https://ratatui.rs/)** - The excellent Rust TUI framework powering this app
- **[Awesome ratatui](https://github.com/ratatui/awesome-ratatui)** - Inspiration from the ecosystem, especially:
  - [termscp](https://github.com/veeso/termscp) - Terminal file transfer
  - [tgv](https://github.com/chroma-pipeline/tgv) - Terminal genome viewer

SeqTUI handles frameshifts consistently with [MACSE](https://github.com/ranwez/MACSE_V2_PIPELINES)

Alternatives to SeqTUI include:
- **[termal](https://github.com/sib-swiss/termal)** - a terminal alignment viewer with conservation scores 
- **[alen](https://github.com/jakobnissen/alen)** - a simpler terminal alignment viewer
- **[alv](https://github.com/arvestad/alv)** - a  terminal alignment viewer to be combined with linux commands
- **[segul](https://www.segul.app)** - a powerful phylogenomic tools (no alignment viewer) 

## License
MIT
