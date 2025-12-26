//! Data model for the alignment viewer.
//!
//! This module contains all data structures for representing:
//! - Sequences and alignments
//! - Viewport state
//! - Application state
//! - Loading state for async operations
//!
//! The design allows for future extensions like filtering, codon views,
//! and translation between nucleotides and amino acids.

use std::ops::Range;
use std::path::PathBuf;

/// Type of biological sequence with nucleotide ratio.
/// The ratio indicates the proportion of nucleotide characters (ACGTUN) found.
/// - ratio > 0.8: displayed with nucleotide colors
/// - ratio < 0.5: likely amino acid, error for NT-requiring operations
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SequenceType {
    /// Ratio of nucleotide characters (0.0 to 1.0)
    pub nt_ratio: f64,
}

impl Default for SequenceType {
    fn default() -> Self {
        Self { nt_ratio: 1.0 }
    }
}

impl SequenceType {
    /// Creates a new SequenceType with the given NT ratio.
    pub fn new(nt_ratio: f64) -> Self {
        Self { nt_ratio }
    }

    /// Returns true if this appears to be nucleotide (>80% NT chars).
    /// Used for display coloring.
    #[inline]
    pub fn is_nucleotide(&self) -> bool {
        self.nt_ratio > 0.8
    }

    /// Returns true if this appears to be amino acid (≤80% NT chars).
    #[inline]
    pub fn is_amino_acid(&self) -> bool {
        !self.is_nucleotide()
    }

    /// Returns true if this is likely NOT nucleotide (<50% NT chars).
    /// Used for error checking in NT-requiring operations (translation, VCF).
    #[inline]
    pub fn is_likely_not_nucleotide(&self) -> bool {
        self.nt_ratio < 0.5
    }

    // Legacy constants for compatibility
    /// Nucleotide sequence type (100% NT ratio)
    pub const NUCLEOTIDE: SequenceType = SequenceType { nt_ratio: 1.0 };
    /// Amino acid sequence type (0% NT ratio)
    pub const AMINO_ACID: SequenceType = SequenceType { nt_ratio: 0.0 };
}

/// Loading state for async operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LoadingState {
    /// No loading in progress, alignment is ready
    Ready,
    /// Loading file from disk
    LoadingFile {
        /// File path being loaded
        path: PathBuf,
        /// Progress message (e.g., "Parsing NEXUS file...")
        message: String,
        /// Number of sequences loaded so far (if available)
        sequences_loaded: Option<usize>,
    },
    /// Translating sequences
    Translating {
        /// Progress message
        message: String,
        /// Number of sequences translated so far
        sequences_done: usize,
        /// Total number of sequences
        total: usize,
    },
}

impl Default for LoadingState {
    fn default() -> Self {
        LoadingState::Ready
    }
}

impl LoadingState {
    /// Returns true if currently loading/processing.
    pub fn is_loading(&self) -> bool {
        !matches!(self, LoadingState::Ready)
    }

    /// Returns the display message for the loading state.
    pub fn message(&self) -> Option<&str> {
        match self {
            LoadingState::Ready => None,
            LoadingState::LoadingFile { message, .. } => Some(message),
            LoadingState::Translating { message, .. } => Some(message),
        }
    }

    /// Returns progress as a fraction (0.0 to 1.0) if available.
    pub fn progress(&self) -> Option<f64> {
        match self {
            LoadingState::Ready => None,
            LoadingState::LoadingFile { .. } => None, // Indeterminate
            LoadingState::Translating { sequences_done, total, .. } => {
                if *total > 0 {
                    Some(*sequences_done as f64 / *total as f64)
                } else {
                    None
                }
            }
        }
    }
}

/// A file entry in the file browser.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileEntry {
    /// File or directory name
    pub name: String,
    /// Full path
    pub path: PathBuf,
    /// Whether this is a directory
    pub is_dir: bool,
}

/// State for the file browser popup.
#[derive(Debug, Clone)]
pub struct FileBrowserState {
    /// Current directory being browsed
    pub current_dir: PathBuf,
    /// List of entries in the current directory
    pub entries: Vec<FileEntry>,
    /// Currently selected index (0-based)
    pub selected: usize,
    /// Scroll offset for long lists
    pub scroll_offset: usize,
    /// Error message that triggered the browser (shown in title)
    pub error_message: String,
}

impl FileBrowserState {
    /// Creates a new file browser starting at the given directory.
    pub fn new(start_dir: PathBuf, error_message: String) -> Self {
        let mut browser = Self {
            current_dir: start_dir,
            entries: Vec::new(),
            selected: 0,
            scroll_offset: 0,
            error_message,
        };
        browser.refresh_entries();
        browser
    }

    /// Refreshes the list of entries from the current directory.
    pub fn refresh_entries(&mut self) {
        self.entries.clear();
        self.selected = 0;
        self.scroll_offset = 0;

        // Add parent directory entry if not at root
        if let Some(parent) = self.current_dir.parent() {
            self.entries.push(FileEntry {
                name: "..".to_string(),
                path: parent.to_path_buf(),
                is_dir: true,
            });
        }

        // Read directory contents
        if let Ok(read_dir) = std::fs::read_dir(&self.current_dir) {
            let mut dirs: Vec<FileEntry> = Vec::new();
            let mut files: Vec<FileEntry> = Vec::new();

            for entry in read_dir.flatten() {
                let path = entry.path();
                let name = entry.file_name().to_string_lossy().to_string();
                
                // Skip hidden files (starting with .)
                if name.starts_with('.') {
                    continue;
                }

                let is_dir = path.is_dir();
                let entry = FileEntry { name, path, is_dir };

                if is_dir {
                    dirs.push(entry);
                } else {
                    // Only show sequence files
                    let ext = entry.path.extension()
                        .and_then(|e| e.to_str())
                        .map(|e| e.to_lowercase());
                    
                    if matches!(ext.as_deref(), 
                        Some("fasta" | "fa" | "fna" | "faa" | "fas" | 
                             "phy" | "phylip" | 
                             "nex" | "nexus" | "nxs")) {
                        files.push(entry);
                    }
                }
            }

            // Sort directories and files separately, then combine
            dirs.sort_by(|a, b| a.name.to_lowercase().cmp(&b.name.to_lowercase()));
            files.sort_by(|a, b| a.name.to_lowercase().cmp(&b.name.to_lowercase()));

            self.entries.extend(dirs);
            self.entries.extend(files);
        }
    }

    /// Moves selection up.
    pub fn select_prev(&mut self) {
        if self.selected > 0 {
            self.selected -= 1;
            // Adjust scroll if needed
            if self.selected < self.scroll_offset {
                self.scroll_offset = self.selected;
            }
        }
    }

    /// Moves selection down.
    pub fn select_next(&mut self) {
        if self.selected + 1 < self.entries.len() {
            self.selected += 1;
        }
    }

    /// Adjusts scroll offset for the given visible height.
    pub fn adjust_scroll(&mut self, visible_height: usize) {
        if self.selected >= self.scroll_offset + visible_height {
            self.scroll_offset = self.selected - visible_height + 1;
        }
    }

    /// Returns the currently selected entry, if any.
    pub fn selected_entry(&self) -> Option<&FileEntry> {
        self.entries.get(self.selected)
    }

    /// Enters the selected directory.
    pub fn enter_selected(&mut self) -> Option<PathBuf> {
        if let Some(entry) = self.selected_entry() {
            if entry.is_dir {
                self.current_dir = entry.path.clone();
                self.refresh_entries();
                None
            } else {
                // Return the selected file path
                Some(entry.path.clone())
            }
        } else {
            None
        }
    }
}

/// Represents a single sequence with its identifier and data.
/// 
/// Sequence data is stored as `Vec<u8>` (ASCII bytes) rather than `String`
/// for efficiency - biological sequences only use ASCII characters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    /// The sequence identifier (from FASTA header, without '>')
    pub id: String,
    /// The sequence data as ASCII bytes (nucleotides or amino acids)
    data: Vec<u8>,
}

impl Sequence {
    /// Creates a new sequence from a string.
    pub fn new(id: impl Into<String>, data: impl AsRef<str>) -> Self {
        Self {
            id: id.into(),
            data: data.as_ref().as_bytes().to_vec(),
        }
    }

    /// Creates a new sequence from raw bytes (zero-copy when possible).
    pub fn from_bytes(id: impl Into<String>, data: Vec<u8>) -> Self {
        Self {
            id: id.into(),
            data,
        }
    }

    /// Returns the length of the sequence.
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Returns true if the sequence is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Gets a character at a specific position (O(1) direct byte access).
    #[inline]
    pub fn char_at(&self, pos: usize) -> Option<char> {
        self.data.get(pos).map(|&b| b as char)
    }

    /// Gets a byte at a specific position.
    #[inline]
    pub fn byte_at(&self, pos: usize) -> Option<u8> {
        self.data.get(pos).copied()
    }

    /// Gets the raw bytes of the sequence.
    #[inline]
    pub fn as_bytes(&self) -> &[u8] {
        &self.data
    }

    /// Takes ownership of the sequence data (for thread transfer).
    /// This avoids copying when the sequence won't be used again.
    #[inline]
    pub fn into_bytes(self) -> Vec<u8> {
        self.data
    }

    /// Clones the sequence data (when ownership transfer isn't possible).
    #[inline]
    pub fn clone_bytes(&self) -> Vec<u8> {
        self.data.clone()
    }

    /// Gets a slice of the sequence data as a string.
    /// This is safe because biological sequences are always ASCII.
    pub fn slice(&self, range: Range<usize>) -> &str {
        let start = range.start.min(self.data.len());
        let end = range.end.min(self.data.len());
        // SAFETY: Biological sequences are always valid ASCII/UTF-8
        unsafe { std::str::from_utf8_unchecked(&self.data[start..end]) }
    }

    /// Gets the entire sequence as a string slice.
    #[inline]
    pub fn as_str(&self) -> &str {
        // SAFETY: Biological sequences are always valid ASCII/UTF-8
        unsafe { std::str::from_utf8_unchecked(&self.data) }
    }
}

/// Represents an alignment of multiple sequences.
#[derive(Debug, Clone)]
pub struct Alignment {
    /// All sequences in the alignment
    pub sequences: Vec<Sequence>,
    /// The expected length of all sequences (if aligned)
    alignment_length: Option<usize>,
    /// Whether all sequences have the same length
    pub is_valid_alignment: bool,
    /// Warning message if sequences have different lengths
    pub warning: Option<String>,
    /// Detected sequence type (nucleotide or amino acid)
    pub sequence_type: SequenceType,
}

impl Alignment {
    /// Creates a new alignment from a vector of sequences.
    pub fn new(sequences: Vec<Sequence>) -> Self {
        let (is_valid, alignment_length, warning) = Self::validate_alignment(&sequences);
        let sequence_type = Self::detect_sequence_type(&sequences);
        Self {
            sequences,
            alignment_length,
            is_valid_alignment: is_valid,
            warning,
            sequence_type,
        }
    }

    /// Detects sequence type by sampling characters across sequences.
    /// 
    /// Samples up to 2500 characters total, taking at most 500 characters per sequence.
    /// This ensures we sample across at least 5 sequences (if available) rather than
    /// just the first sequence, providing a more reliable estimate.
    /// 
    /// Returns a SequenceType with the computed NT ratio.
    fn detect_sequence_type(sequences: &[Sequence]) -> SequenceType {
        if sequences.is_empty() {
            return SequenceType::NUCLEOTIDE;
        }

        let mut nucleotide_chars = 0usize;
        let mut total_chars = 0usize;
        const MAX_CHARS_PER_SEQ: usize = 500;
        const MAX_TOTAL_CHARS: usize = 2500;

        for seq in sequences {
            let mut seq_chars = 0usize;
            for &b in seq.as_bytes() {
                let upper = b.to_ascii_uppercase();
                // Skip gaps, missing data, and spaces
                if matches!(upper, b'-' | b'.' | b' ' | b'?' | b'N') {
                    continue;
                }
                total_chars += 1;
                seq_chars += 1;
                // Nucleotide characters (ACGTU)
                if matches!(upper, b'A' | b'C' | b'G' | b'T' | b'U') {
                    nucleotide_chars += 1;
                }
                // Stop sampling this sequence after MAX_CHARS_PER_SEQ
                if seq_chars >= MAX_CHARS_PER_SEQ {
                    break;
                }
            }
            // Stop sampling after MAX_TOTAL_CHARS total
            if total_chars >= MAX_TOTAL_CHARS {
                break;
            }
        }

        if total_chars == 0 {
            return SequenceType::NUCLEOTIDE;
        }

        SequenceType::new(nucleotide_chars as f64 / total_chars as f64)
    }

    /// Validates that all sequences have the same length.
    fn validate_alignment(sequences: &[Sequence]) -> (bool, Option<usize>, Option<String>) {
        if sequences.is_empty() {
            return (true, None, None);
        }

        let first_len = sequences[0].len();
        let all_same = sequences.iter().all(|s| s.len() == first_len);

        if all_same {
            (true, Some(first_len), None)
        } else {
            let min_len = sequences.iter().map(|s| s.len()).min().unwrap_or(0);
            let max_len = sequences.iter().map(|s| s.len()).max().unwrap_or(0);
            let warning = format!(
                "Warning: Sequences have different lengths (min: {}, max: {}). Not a valid alignment.",
                min_len, max_len
            );
            (false, Some(max_len), Some(warning))
        }
    }

    /// Returns the number of sequences.
    pub fn sequence_count(&self) -> usize {
        self.sequences.len()
    }

    /// Returns the alignment length (max sequence length).
    pub fn alignment_length(&self) -> usize {
        self.alignment_length.unwrap_or(0)
    }

    /// Returns the maximum identifier length (for display purposes).
    pub fn max_id_length(&self) -> usize {
        self.sequences.iter().map(|s| s.id.len()).max().unwrap_or(0)
    }

    /// Gets a sequence by index.
    pub fn get(&self, index: usize) -> Option<&Sequence> {
        self.sequences.get(index)
    }

    /// Returns true if the alignment is empty.
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }
}

/// The viewport defines what portion of the alignment is currently visible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Viewport {
    /// Index of the first visible sequence (row)
    pub first_row: usize,
    /// Index of the first visible column
    pub first_col: usize,
    /// Number of visible rows
    pub visible_rows: usize,
    /// Number of visible columns
    pub visible_cols: usize,
}

impl Viewport {
    /// Creates a new viewport.
    pub fn new(visible_rows: usize, visible_cols: usize) -> Self {
        Self {
            first_row: 0,
            first_col: 0,
            visible_rows,
            visible_cols,
        }
    }

    /// Updates the viewport dimensions.
    pub fn resize(&mut self, visible_rows: usize, visible_cols: usize) {
        self.visible_rows = visible_rows;
        self.visible_cols = visible_cols;
    }

    /// Returns the range of visible rows.
    pub fn row_range(&self) -> Range<usize> {
        self.first_row..self.first_row + self.visible_rows
    }

    /// Returns the range of visible columns.
    pub fn col_range(&self) -> Range<usize> {
        self.first_col..self.first_col + self.visible_cols
    }

    /// Checks if a column is visible.
    pub fn is_col_visible(&self, col: usize) -> bool {
        col >= self.first_col && col < self.first_col + self.visible_cols
    }

    /// Checks if a row is visible.
    pub fn is_row_visible(&self, row: usize) -> bool {
        row >= self.first_row && row < self.first_row + self.visible_rows
    }
}

/// The current cursor position in the alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Cursor {
    /// Current row (sequence index)
    pub row: usize,
    /// Current column (position in sequence)
    pub col: usize,
}

impl Cursor {
    /// Creates a new cursor at origin.
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates a cursor at a specific position.
    pub fn at(row: usize, col: usize) -> Self {
        Self { row, col }
    }
}

/// Help tab sections
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HelpTab {
    #[default]
    Basics,
    Navigation,
    VimNav,
    Search,
    Translation,
}

impl HelpTab {
    /// Returns the next tab (wrapping around)
    pub fn next(self) -> Self {
        match self {
            HelpTab::Basics => HelpTab::Navigation,
            HelpTab::Navigation => HelpTab::VimNav,
            HelpTab::VimNav => HelpTab::Search,
            HelpTab::Search => HelpTab::Translation,
            HelpTab::Translation => HelpTab::Basics,
        }
    }

    /// Returns the previous tab (wrapping around)
    pub fn prev(self) -> Self {
        match self {
            HelpTab::Basics => HelpTab::Translation,
            HelpTab::Navigation => HelpTab::Basics,
            HelpTab::VimNav => HelpTab::Navigation,
            HelpTab::Search => HelpTab::VimNav,
            HelpTab::Translation => HelpTab::Search,
        }
    }

    /// Returns the tab name for display
    pub fn name(self) -> &'static str {
        match self {
            HelpTab::Basics => "Basics",
            HelpTab::Navigation => "Arrow Nav",
            HelpTab::VimNav => "Vim Nav",
            HelpTab::Search => "Search",
            HelpTab::Translation => "Translation",
        }
    }

    /// Returns all tabs in order
    pub fn all() -> &'static [HelpTab] {
        &[HelpTab::Basics, HelpTab::Navigation, HelpTab::VimNav, HelpTab::Search, HelpTab::Translation]
    }
}

/// Application mode for handling different input states.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum AppMode {
    /// Normal navigation mode
    #[default]
    Normal,
    /// Command input mode (after pressing ':')
    Command(String),
    /// Search forward mode (after pressing '/')
    Search(String),
    /// Search backward mode (after pressing '?')
    SearchBackward(String),
    /// Translation settings mode (selecting genetic code and frame)
    TranslationSettings,
}

/// View mode for the alignment (nucleotide or translated amino acid).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ViewMode {
    /// Original nucleotide view
    #[default]
    Nucleotide,
    /// Translated amino acid view
    AminoAcid,
}

/// Translation settings for NT to AA conversion.
#[derive(Debug, Clone)]
pub struct TranslationSettings {
    /// Selected genetic code ID (1-33)
    pub genetic_code_id: u8,
    /// Reading frame (0, 1, or 2 for +1, +2, +3)
    pub frame: usize,
    /// Currently selected genetic code index in the list (for UI)
    pub selected_code_index: usize,
    /// Currently selected frame in the UI (0, 1, 2)
    pub selected_frame: usize,
    /// Scroll offset for the genetic code list
    pub scroll_offset: usize,
    /// Whether user has translated at least once (to show settings on first use)
    pub has_translated: bool,
}

impl Default for TranslationSettings {
    fn default() -> Self {
        Self {
            genetic_code_id: 1, // Standard code
            frame: 0,          // +1 frame
            selected_code_index: 0,
            selected_frame: 0,
            scroll_offset: 0,
            has_translated: false,
        }
    }
}

/// The complete application state.
#[derive(Debug)]
pub struct AppState {
    /// File name (basename without extension) for display
    pub file_name: String,
    /// The original loaded alignment (always nucleotide if translated)
    pub alignment: Alignment,
    /// The translated alignment (if viewing as AA) - lazily computed and cached
    translated_alignment: Option<Alignment>,
    /// Genetic code ID used for the cached translation (to detect if recomputation needed)
    cached_translation_code_id: Option<u8>,
    /// Reading frame used for the cached translation (to detect if recomputation needed)
    cached_translation_frame: Option<usize>,
    /// Current view mode
    pub view_mode: ViewMode,
    /// Translation settings
    pub translation_settings: TranslationSettings,
    /// Current viewport
    pub viewport: Viewport,
    /// Current cursor position
    pub cursor: Cursor,
    /// Current application mode
    pub mode: AppMode,
    /// Whether the application should quit
    pub should_quit: bool,
    /// Status message to display
    pub status_message: Option<String>,
    /// Last search pattern
    pub last_search: Option<String>,
    /// Last search direction (true = backward, false = forward)
    pub last_search_backward: bool,
    /// Whether to show the help overlay
    pub show_help: bool,
    /// Current help tab
    pub help_tab: HelpTab,
    /// Pending 'g' key for g-commands
    pub pending_g: bool,
    /// Pending 'z' key for z-commands
    pub pending_z: bool,
    /// Number buffer for <number>| command
    pub number_buffer: String,
    /// Loading state for async operations
    pub loading_state: LoadingState,
    /// Spinner animation frame (0-3)
    pub spinner_frame: usize,
    /// Error popup to display (shown as centered popup, dismissable with any key)
    pub error_popup: Option<String>,
    /// File browser state (shown when file not found)
    pub file_browser: Option<FileBrowserState>,
    /// Whether to use fancy UI glyphs (non-Windows only)
    pub fancy_ui: bool,
}

impl AppState {
        /// Goes to a specific row (1-indexed for user, like Vim :line command).
        /// Clamps column to alignment length. Returns true if successful, false if out of range.
        fn goto_row_1indexed(&mut self, row: usize) -> bool {
            let num_seqs = self.active_alignment().sequence_count();
            if row > 0 && row <= num_seqs {
                self.cursor.row = row - 1; // 1-indexed for user
                // Keep cursor column, but clamp to alignment length
                let aln_len = self.active_alignment().alignment_length();
                self.cursor.col = self.cursor.col.min(aln_len.saturating_sub(1));
                self.ensure_cursor_visible();
                self.status_message = None; // Clear any previous error
                true
            } else {
                self.status_message = Some(format!("Invalid sequence number: {}", row));
                false
            }
        }
    /// Creates a new application state with the given alignment.
    pub fn new(alignment: Alignment, file_name: String) -> Self {
        let warning = alignment.warning.clone();
        Self {
            file_name,
            alignment,
            translated_alignment: None,
            cached_translation_code_id: None,
            cached_translation_frame: None,
            view_mode: ViewMode::Nucleotide,
            translation_settings: TranslationSettings::default(),
            viewport: Viewport::new(0, 0),
            cursor: Cursor::new(),
            mode: AppMode::Normal,
            should_quit: false,
            status_message: warning,
            last_search: None,
            last_search_backward: false,
            show_help: false,
            help_tab: HelpTab::default(),
            pending_g: false,
            pending_z: false,
            number_buffer: String::new(),
            loading_state: LoadingState::Ready,
            spinner_frame: 0,
            error_popup: None,
            file_browser: None,
            fancy_ui: false,
        }
    }

    /// Creates an empty application state for async loading.
    pub fn new_loading(file_name: String, path: PathBuf) -> Self {
        Self {
            file_name: file_name.clone(),
            alignment: Alignment::new(vec![]),
            translated_alignment: None,
            cached_translation_code_id: None,
            cached_translation_frame: None,
            view_mode: ViewMode::Nucleotide,
            translation_settings: TranslationSettings::default(),
            viewport: Viewport::new(0, 0),
            cursor: Cursor::new(),
            mode: AppMode::Normal,
            should_quit: false,
            status_message: None,
            last_search: None,
            last_search_backward: false,
            show_help: false,
            help_tab: HelpTab::default(),
            pending_g: false,
            pending_z: false,
            number_buffer: String::new(),
            loading_state: LoadingState::LoadingFile {
                path,
                message: format!("Loading {}...", file_name),
                sequences_loaded: None,
            },
            spinner_frame: 0,
            error_popup: None,
            file_browser: None,
            fancy_ui: false,
        }
    }

    /// Updates the alignment after async loading completes.
    pub fn set_alignment(&mut self, alignment: Alignment) {
        let warning = alignment.warning.clone();
        self.alignment = alignment;
        self.loading_state = LoadingState::Ready;
        if let Some(w) = warning {
            self.status_message = Some(w);
        }
    }

    /// Sets an error state after loading fails.
    /// If the error is "file not found", opens the file browser.
    pub fn set_loading_error(&mut self, error: String, file_path: Option<PathBuf>) {
        self.loading_state = LoadingState::Ready;
        
        // Check if this is a "file not found" error
        let is_not_found = error.contains("No such file") || error.contains("not found");
        
        if is_not_found {
            // Open file browser starting from the file's directory or current directory
            let start_dir = file_path
                .as_ref()
                .and_then(|p| p.parent())
                .filter(|p| !p.as_os_str().is_empty()) // Filter out empty parent paths
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")));
            
            let file_name = file_path
                .as_ref()
                .and_then(|p| p.file_name())
                .and_then(|n| n.to_str())
                .unwrap_or("file");
            
            let error_msg = format!("File not found: {}", file_name);
            self.file_browser = Some(FileBrowserState::new(start_dir, error_msg));
        } else {
            // Show error popup for other errors
            self.error_popup = Some(error);
        }
    }

    /// Shows an error popup with the given message.
    pub fn show_error_popup(&mut self, message: String) {
        self.error_popup = Some(message);
    }

    /// Dismisses the error popup.
    pub fn dismiss_error_popup(&mut self) {
        self.error_popup = None;
    }

    /// Closes the file browser.
    pub fn close_file_browser(&mut self) {
        self.file_browser = None;
    }

    /// Moves file browser selection up.
    pub fn file_browser_up(&mut self) {
        if let Some(browser) = &mut self.file_browser {
            browser.select_prev();
        }
    }

    /// Moves file browser selection down.
    pub fn file_browser_down(&mut self) {
        if let Some(browser) = &mut self.file_browser {
            browser.select_next();
        }
    }

    /// Selects the current file browser entry.
    /// Returns Some(path) if a file was selected, None if navigating into directory.
    pub fn file_browser_select(&mut self) -> Option<PathBuf> {
        if let Some(browser) = &mut self.file_browser {
            browser.enter_selected()
        } else {
            None
        }
    }

    /// Goes to parent directory in file browser.
    pub fn file_browser_parent(&mut self) {
        if let Some(browser) = &mut self.file_browser {
            if let Some(parent) = browser.current_dir.parent() {
                browser.current_dir = parent.to_path_buf();
                browser.refresh_entries();
            }
        }
    }

    /// Quits the file browser. If an alignment is loaded, just close the browser.
    /// If no alignment is loaded, quit the application.
    pub fn file_browser_quit(&mut self) {
        self.file_browser = None;
        // Only quit if no alignment is loaded
        if self.alignment.sequence_count() == 0 {
            self.should_quit = true;
        }
    }

    /// Sets the translated alignment after async translation completes.
    pub fn set_translated_alignment(&mut self, alignment: Alignment) {
        self.translated_alignment = Some(alignment);
        // Cache the settings used for this translation
        self.cached_translation_code_id = Some(self.translation_settings.genetic_code_id);
        self.cached_translation_frame = Some(self.translation_settings.frame);
        self.loading_state = LoadingState::Ready;
        self.view_mode = ViewMode::AminoAcid;
        
        // Convert NT cursor position to AA position
        let frame = self.translation_settings.frame;
        let aa_col = if self.cursor.col >= frame {
            (self.cursor.col - frame) / 3
        } else {
            0
        };
        
        let aa_len = self.translated_alignment.as_ref()
            .map(|a| a.alignment_length())
            .unwrap_or(0);
        self.cursor.col = aa_col.min(aa_len.saturating_sub(1));
        
        self.ensure_cursor_visible();
        self.status_message = Some(format!(
            "Translated using code {} (frame +{})",
            self.translation_settings.genetic_code_id,
            self.translation_settings.frame + 1
        ));
    }

    /// Returns true if we should start background translation.
    /// Called when user triggers translation command.
    /// Returns false if we already have a cached translation with matching settings.
    pub fn should_start_translation(&self) -> bool {
        self.alignment.sequence_type.is_nucleotide()
            && !self.loading_state.is_loading()
            && !self.has_valid_cached_translation()
    }
    
    /// Returns true if we have a cached translated alignment that matches
    /// the current translation settings (code_id and frame).
    pub fn has_valid_cached_translation(&self) -> bool {
        if self.translated_alignment.is_none() {
            return false;
        }
        let code_matches = self.cached_translation_code_id == Some(self.translation_settings.genetic_code_id);
        let frame_matches = self.cached_translation_frame == Some(self.translation_settings.frame);
        code_matches && frame_matches
    }
    
    /// Switches to the cached amino acid view without recomputation.
    /// Call this when has_valid_cached_translation() returns true.
    pub fn switch_to_cached_aa_view(&mut self) {
        if self.translated_alignment.is_none() {
            return;
        }
        self.view_mode = ViewMode::AminoAcid;
        
        // Convert NT cursor position to AA position
        let frame = self.translation_settings.frame;
        let aa_col = if self.cursor.col >= frame {
            (self.cursor.col - frame) / 3
        } else {
            0
        };
        
        let aa_len = self.translated_alignment.as_ref()
            .map(|a| a.alignment_length())
            .unwrap_or(0);
        self.cursor.col = aa_col.min(aa_len.saturating_sub(1));
        
        self.ensure_cursor_visible();
        self.status_message = Some(format!(
            "Using cached translation (code {}, frame +{})",
            self.cached_translation_code_id.unwrap_or(1),
            self.cached_translation_frame.unwrap_or(0) + 1
        ));
    }

    /// Advances the spinner animation frame.
    pub fn tick_spinner(&mut self) {
        self.spinner_frame = (self.spinner_frame + 1) % 4;
    }

    /// Returns the current spinner character.
    pub fn spinner_char(&self) -> char {
        const SPINNER: [char; 4] = ['|', '/', '-', '\\'];
        SPINNER[self.spinner_frame]
    }

    /// Returns the currently active alignment (original or translated).
    pub fn active_alignment(&self) -> &Alignment {
        match self.view_mode {
            ViewMode::Nucleotide => &self.alignment,
            ViewMode::AminoAcid => self.translated_alignment.as_ref().unwrap_or(&self.alignment),
        }
    }

    /// Returns whether a translated alignment exists.
    pub fn has_translated_alignment(&self) -> bool {
        self.translated_alignment.is_some()
    }

    /// Updates the viewport size based on terminal dimensions.
    pub fn update_viewport_size(&mut self, rows: usize, cols: usize) {
        self.viewport.resize(rows, cols);
        self.ensure_cursor_visible();
    }

    /// Moves the cursor up by one row.
    pub fn move_up(&mut self) {
        if self.cursor.row > 0 {
            self.cursor.row -= 1;
            self.ensure_cursor_visible();
        }
    }

    /// Moves the cursor down by one row.
    pub fn move_down(&mut self) {
        if self.cursor.row + 1 < self.active_alignment().sequence_count() {
            self.cursor.row += 1;
            self.ensure_cursor_visible();
        }
    }

    /// Moves the cursor left by one column.
    pub fn move_left(&mut self) {
        if self.cursor.col > 0 {
            self.cursor.col -= 1;
            self.ensure_cursor_visible();
        }
    }

    /// Moves the cursor right by one column.
    pub fn move_right(&mut self) {
        if self.cursor.col + 1 < self.active_alignment().alignment_length() {
            self.cursor.col += 1;
            self.ensure_cursor_visible();
        }
    }

    /// Moves the cursor up by half a page.
    pub fn half_page_up(&mut self) {
        let half_page = self.viewport.visible_rows / 2;
        self.cursor.row = self.cursor.row.saturating_sub(half_page);
        self.ensure_cursor_visible();
    }

    /// Moves the cursor down by half a page.
    pub fn half_page_down(&mut self) {
        let half_page = self.viewport.visible_rows / 2;
        let max_row = self.active_alignment().sequence_count().saturating_sub(1);
        self.cursor.row = (self.cursor.row + half_page).min(max_row);
        self.ensure_cursor_visible();
    }

    /// Moves the cursor left by half a screen width.
    pub fn half_page_left(&mut self) {
        self.clear_pending();
        let half_page = self.viewport.visible_cols / 2;
        self.cursor.col = self.cursor.col.saturating_sub(half_page);
        self.ensure_cursor_visible();
    }

    /// Moves the cursor right by half a screen width.
    pub fn half_page_right(&mut self) {
        self.clear_pending();
        let half_page = self.viewport.visible_cols / 2;
        let max_col = self.active_alignment().alignment_length().saturating_sub(1);
        self.cursor.col = (self.cursor.col + half_page).min(max_col);
        self.ensure_cursor_visible();
    }

    /// Moves the cursor up by a full page.
    pub fn page_up(&mut self) {
        let page = self.viewport.visible_rows.saturating_sub(1).max(1);
        self.cursor.row = self.cursor.row.saturating_sub(page);
        self.ensure_cursor_visible();
    }

    /// Moves the cursor down by a full page.
    pub fn page_down(&mut self) {
        let page = self.viewport.visible_rows.saturating_sub(1).max(1);
        let max_row = self.active_alignment().sequence_count().saturating_sub(1);
        self.cursor.row = (self.cursor.row + page).min(max_row);
        self.ensure_cursor_visible();
    }

    /// Returns true if a character is a word delimiter (gap or stop codon).
    #[inline]
    fn is_word_delimiter(c: u8) -> bool {
        c == b'-' || c == b'*'
    }

    /// Moves to the start of the next word (w).
    /// Words are delimited by gaps (-) and stops (*).
    pub fn word_forward(&mut self) {
        let alignment = self.active_alignment();
        let aln_len = alignment.alignment_length();
        if aln_len == 0 || self.cursor.row >= alignment.sequence_count() {
            return;
        }

        let seq = alignment.sequences[self.cursor.row].as_bytes();
        let mut col = self.cursor.col;

        // If on a non-delimiter, skip to end of current word
        if col < seq.len() && !Self::is_word_delimiter(seq[col]) {
            while col < seq.len() && !Self::is_word_delimiter(seq[col]) {
                col += 1;
            }
        }

        // Skip delimiters to find start of next word
        while col < seq.len() && Self::is_word_delimiter(seq[col]) {
            col += 1;
        }

        // Clamp to valid range
        self.cursor.col = col.min(aln_len.saturating_sub(1));
        self.ensure_cursor_visible();
    }

    /// Moves to the start of the previous word (b).
    /// Words are delimited by gaps (-) and stops (*).
    pub fn word_backward(&mut self) {
        let alignment = self.active_alignment();
        if alignment.alignment_length() == 0 || self.cursor.row >= alignment.sequence_count() {
            return;
        }

        let seq = alignment.sequences[self.cursor.row].as_bytes();
        let mut col = self.cursor.col;

        // If at start, nothing to do
        if col == 0 {
            return;
        }

        // Move back one position first
        col -= 1;

        // Skip delimiters backward
        while col > 0 && Self::is_word_delimiter(seq[col]) {
            col -= 1;
        }

        // Now we're on a word character (or at position 0)
        // Find the start of this word
        while col > 0 && !Self::is_word_delimiter(seq[col - 1]) {
            col -= 1;
        }

        self.cursor.col = col;
        self.ensure_cursor_visible();
    }

    /// Moves to the end of the current/next word (e).
    /// Words are delimited by gaps (-) and stops (*).
    pub fn word_end(&mut self) {
        let alignment = self.active_alignment();
        let aln_len = alignment.alignment_length();
        if aln_len == 0 || self.cursor.row >= alignment.sequence_count() {
            return;
        }

        let seq = alignment.sequences[self.cursor.row].as_bytes();
        let mut col = self.cursor.col;

        // Move forward at least one position
        if col + 1 < seq.len() {
            col += 1;
        }

        // Skip delimiters
        while col < seq.len() && Self::is_word_delimiter(seq[col]) {
            col += 1;
        }

        // Find end of this word
        while col + 1 < seq.len() && !Self::is_word_delimiter(seq[col + 1]) {
            col += 1;
        }

        // Clamp to valid range
        self.cursor.col = col.min(aln_len.saturating_sub(1));
        self.ensure_cursor_visible();
    }

    /// Ensures the cursor is visible in the viewport, with centering behavior.
    fn ensure_cursor_visible(&mut self) {
        // Vertical scrolling - keep cursor in view
        if self.cursor.row < self.viewport.first_row {
            self.viewport.first_row = self.cursor.row;
        } else if self.cursor.row >= self.viewport.first_row + self.viewport.visible_rows {
            self.viewport.first_row = self.cursor.row.saturating_sub(self.viewport.visible_rows - 1);
        }

        // Horizontal scrolling - center when reaching edge
        if self.cursor.col < self.viewport.first_col {
            // Cursor went left of viewport - center it
            self.center_column();
        } else if self.cursor.col >= self.viewport.first_col + self.viewport.visible_cols {
            // Cursor went right of viewport - center it
            self.center_column();
        }

        // Clamp viewport to valid bounds
        self.clamp_viewport();
    }

    /// Centers the current column in the viewport.
    fn center_column(&mut self) {
        if self.viewport.visible_cols > 0 {
            let half = self.viewport.visible_cols / 2;
            self.viewport.first_col = self.cursor.col.saturating_sub(half);
        }
    }

    /// Clamps the viewport to valid alignment bounds.
    fn clamp_viewport(&mut self) {
        let seq_count = self.active_alignment().sequence_count();
        let aln_length = self.active_alignment().alignment_length();
        let max_row = seq_count.saturating_sub(1);
        let max_col = aln_length.saturating_sub(1);

        // Ensure we don't scroll past the end
        if self.viewport.first_row + self.viewport.visible_rows > seq_count {
            self.viewport.first_row = seq_count.saturating_sub(self.viewport.visible_rows);
        }

        if self.viewport.first_col + self.viewport.visible_cols > aln_length {
            self.viewport.first_col = aln_length.saturating_sub(self.viewport.visible_cols);
        }

        // Clamp cursor to valid bounds
        self.cursor.row = self.cursor.row.min(max_row);
        self.cursor.col = self.cursor.col.min(max_col);
    }

    /// Enters command mode.
    pub fn enter_command_mode(&mut self) {
        self.mode = AppMode::Command(String::new());
    }

    /// Handles a character input in command mode.
    pub fn command_input(&mut self, c: char) {
        if let AppMode::Command(ref mut cmd) = self.mode {
            cmd.push(c);
        }
    }

    /// Handles backspace in command mode.
    pub fn command_backspace(&mut self) {
        if let AppMode::Command(ref mut cmd) = self.mode {
            cmd.pop();
            if cmd.is_empty() {
                self.mode = AppMode::Normal;
            }
        }
    }

    /// Executes the current command.
    /// Returns true if translation should be started (handled by controller).
    pub fn execute_command(&mut self) -> bool {
        let mut start_translation = false;
        if let AppMode::Command(ref cmd) = self.mode.clone() {
            match cmd.as_str() {
                "$" => {
                    let num_seqs = self.active_alignment().sequence_count();
                    if num_seqs > 0 {
                        self.goto_row_1indexed(num_seqs);
                    } else {
                        self.status_message = Some("No sequences loaded".to_string());
                    }
                }
                "q" | "quit" => self.should_quit = true,
                "h" | "help" => self.show_help(),
                "e" | "edit" => {
                    // Open file browser to select a new file
                    let start_dir = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
                    self.file_browser = Some(FileBrowserState::new(start_dir, "Open file".to_string()));
                    self.mode = AppMode::Normal;
                    return false;
                }
                "asAA" | "asaa" => {
                    if !self.alignment.sequence_type.is_nucleotide() {
                        self.status_message = Some("Cannot translate: not a nucleotide sequence".to_string());
                    } else if self.has_valid_cached_translation() {
                        // We have a valid cached translation - switch to it directly
                        self.switch_to_cached_aa_view();
                    } else if !self.translation_settings.has_translated {
                        // First time: show settings so user knows the options
                        self.enter_translation_settings();
                        return false; // Don't reset to Normal
                    } else {
                        // Subsequent times: request translation (controller will start background thread)
                        start_translation = true;
                    }
                }
                "asNT" | "asnt" => self.switch_to_nucleotide_view(),
                "setcode" => {
                    self.enter_translation_settings();
                    return false; // Don't reset to Normal - enter_translation_settings sets the mode
                }
                "w" => {
                    // :w alone without filename
                    self.status_message = Some("Usage: :w filename.fasta".to_string());
                }
                _ => {
                    // Handle :w filename - save to FASTA
                    if cmd.starts_with("w ") {
                        let filename = cmd[2..].trim();
                        if filename.is_empty() {
                            self.status_message = Some("Usage: :w filename.fasta".to_string());
                        } else {
                            match self.write_fasta(filename) {
                                Ok(count) => {
                                    self.status_message = Some(format!(
                                        "Saved {} sequences to {}", count, filename
                                    ));
                                }
                                Err(e) => {
                                    self.status_message = Some(format!("Error saving: {}", e));
                                }
                            }
                        }
                    }
                    // Handle :number for row/sequence navigation (like Vim's :line)
                    else if let Ok(row) = cmd.parse::<usize>() {
                        let num_seqs = self.active_alignment().sequence_count();
                        if num_seqs > 0 {
                            self.goto_row_1indexed(row);
                        } else {
                            self.status_message = Some("No sequences loaded".to_string());
                        }
                    } else {
                        self.status_message = Some(format!("Unknown command: {}", cmd));
                    }
                }
            }
        }
        self.mode = AppMode::Normal;
        start_translation
    }

    /// Enters translation settings mode.
    pub fn enter_translation_settings(&mut self) {
        // Only allow translation if the original sequence is nucleotide
        if !self.alignment.sequence_type.is_nucleotide() {
            self.status_message = Some("Cannot translate: not a nucleotide sequence".to_string());
            return;
        }
        self.mode = AppMode::TranslationSettings;
    }

    /// Writes the current view (NT or AA) to a FASTA file.
    /// Sequences are written on a single line (convenient for bash processing).
    /// Returns the number of sequences written.
    pub fn write_fasta(&self, filename: &str) -> std::io::Result<usize> {
        use std::io::Write;
        
        let alignment = self.active_alignment();
        let mut file = std::fs::File::create(filename)?;
        
        for seq in &alignment.sequences {
            // Write header
            writeln!(file, ">{}", seq.id)?;
            // Write sequence on single line
            writeln!(file, "{}", seq.as_str())?;
        }
        
        Ok(alignment.sequence_count())
    }

    /// Switches to nucleotide view.
    pub fn switch_to_nucleotide_view(&mut self) {
        if self.view_mode == ViewMode::Nucleotide {
            self.status_message = Some("Already in nucleotide view".to_string());
            return;
        }
        
        // Convert AA position to NT position: first nucleotide of the codon
        // AA position n corresponds to NT positions (n*3 + frame) to (n*3 + frame + 2)
        let frame = self.translation_settings.frame as usize;
        let nt_col = self.cursor.col * 3 + frame;
        self.cursor.col = nt_col.min(self.alignment.alignment_length().saturating_sub(1));
        
        self.view_mode = ViewMode::Nucleotide;
        
        // Keep the translated alignment cached for fast switching back
        // It will be invalidated if genetic code or frame changes
        
        self.ensure_cursor_visible();
        self.status_message = Some("Switched to nucleotide view".to_string());
    }

    /// Switches to amino acid view (performs translation).
    pub fn switch_to_amino_acid_view(&mut self) {
        use crate::genetic_code::GeneticCodes;
        
        // Perform translation
        let codes = GeneticCodes::new();
        let code = codes.get(self.translation_settings.genetic_code_id)
            .unwrap_or_else(|| codes.default_code());
        let frame = self.translation_settings.frame as usize;
        
        // Translate all sequences (single-threaded - fast enough now)
        let translated_seqs: Vec<crate::model::Sequence> = self.alignment.sequences
            .iter()
            .map(|seq| {
                let aa_data = code.translate_sequence(seq.as_bytes(), frame);
                Sequence::from_bytes(seq.id.clone(), aa_data)
            })
            .collect();
        
        let mut translated = Alignment::new(translated_seqs);
        // Force sequence type to AminoAcid
        translated.sequence_type = SequenceType::AMINO_ACID;
        self.translated_alignment = Some(translated);
        
        // Convert NT cursor position to AA position
        // NT position n in frame f corresponds to AA position (n - f) / 3
        let aa_col = if self.cursor.col >= frame {
            (self.cursor.col - frame) / 3
        } else {
            0
        };
        
        let aa_len = self.translated_alignment.as_ref()
            .map(|a| a.alignment_length())
            .unwrap_or(0);
        self.cursor.col = aa_col.min(aa_len.saturating_sub(1));
        
        self.view_mode = ViewMode::AminoAcid;
        self.ensure_cursor_visible();
        self.status_message = Some(format!(
            "Translated using code {} (frame +{})",
            self.translation_settings.genetic_code_id,
            self.translation_settings.frame + 1
        ));
    }

    /// Confirms translation settings (does NOT start translation - that's handled by controller).
    /// Returns true if translation should be started.
    pub fn confirm_translation_settings(&mut self) -> bool {
        // Copy selected values to actual settings
        self.translation_settings.genetic_code_id = {
            use crate::genetic_code::GeneticCodes;
            let codes = GeneticCodes::new();
            codes.all().get(self.translation_settings.selected_code_index)
                .map(|c| c.id)
                .unwrap_or(1)
        };
        self.translation_settings.frame = self.translation_settings.selected_frame;
        self.translation_settings.has_translated = true;
        
        self.mode = AppMode::Normal;
        
        // Check if we can use cached translation (settings match)
        if self.has_valid_cached_translation() {
            self.switch_to_cached_aa_view();
            return false;
        }
        
        // Return true to signal controller should start background translation
        true
    }

    /// Cancels translation settings.
    pub fn cancel_translation_settings(&mut self) {
        self.mode = AppMode::Normal;
    }

    /// Moves selection up in translation settings.
    pub fn translation_settings_up(&mut self) {
        if self.translation_settings.selected_code_index > 0 {
            self.translation_settings.selected_code_index -= 1;
            // Adjust scroll if needed
            if self.translation_settings.selected_code_index < self.translation_settings.scroll_offset {
                self.translation_settings.scroll_offset = self.translation_settings.selected_code_index;
            }
        }
    }

    /// Moves selection down in translation settings.
    pub fn translation_settings_down(&mut self) {
        use crate::genetic_code::GeneticCodes;
        let codes = GeneticCodes::new();
        let max_index = codes.all().len().saturating_sub(1);
        if self.translation_settings.selected_code_index < max_index {
            self.translation_settings.selected_code_index += 1;
        }
    }

    /// Cycles frame selection left in translation settings.
    pub fn translation_settings_frame_left(&mut self) {
        if self.translation_settings.selected_frame > 0 {
            self.translation_settings.selected_frame -= 1;
        }
    }

    /// Cycles frame selection right in translation settings.
    pub fn translation_settings_frame_right(&mut self) {
        if self.translation_settings.selected_frame < 2 {
            self.translation_settings.selected_frame += 1;
        }
    }

    /// Shows help information.
    fn show_help(&mut self) {
        self.show_help = !self.show_help;
    }

    /// Toggles help off (for any key press while help is shown).
    pub fn dismiss_help(&mut self) {
        self.show_help = false;
    }

    /// Navigates to the next help tab.
    pub fn help_next_tab(&mut self) {
        self.help_tab = self.help_tab.next();
    }

    /// Navigates to the previous help tab.
    pub fn help_prev_tab(&mut self) {
        self.help_tab = self.help_tab.prev();
    }

    /// Clears any pending key state.
    fn clear_pending(&mut self) {
        self.pending_g = false;
        self.pending_z = false;
        self.number_buffer.clear();
    }

    /// Sets the pending 'g' state for g-commands.
    pub fn set_pending_g(&mut self) {
        self.pending_g = true;
        self.pending_z = false;
        self.number_buffer.clear();
    }

    /// Sets the pending 'z' state for z-commands.
    pub fn set_pending_z(&mut self) {
        self.pending_z = true;
        self.pending_g = false;
        self.number_buffer.clear();
    }

    /// Accumulates a digit for the number prefix.
    pub fn accumulate_digit(&mut self, c: char) {
        self.pending_g = false;
        self.pending_z = false;
        self.number_buffer.push(c);
    }

    /// Executes goto column with the accumulated number.
    pub fn execute_goto_column(&mut self) {
        if let Ok(col) = self.number_buffer.parse::<usize>() {
            self.goto_column(col);
        }
        self.clear_pending();
    }

    /// Goes to the first column (0).
    pub fn goto_first_column(&mut self) {
        self.clear_pending();
        self.cursor.col = 0;
        self.ensure_cursor_visible();
    }

    /// Goes to the last column (end of alignment).
    pub fn goto_last_column(&mut self) {
        self.clear_pending();
        let aln_len = self.active_alignment().alignment_length();
        if aln_len > 0 {
            self.cursor.col = aln_len - 1;
            self.ensure_cursor_visible();
        }
    }

    /// Goes to the first visible column (g0).
    pub fn goto_first_visible_column(&mut self) {
        self.clear_pending();
        self.cursor.col = self.viewport.first_col;
    }

    /// Goes to the middle of the visible area (gm).
    pub fn goto_middle_visible_column(&mut self) {
        self.clear_pending();
        let middle = self.viewport.first_col + self.viewport.visible_cols / 2;
        self.cursor.col = middle.min(self.active_alignment().alignment_length().saturating_sub(1));
    }

    /// Goes to the last visible column (g$).
    pub fn goto_last_visible_column(&mut self) {
        self.clear_pending();
        let last_visible = self.viewport.first_col + self.viewport.visible_cols.saturating_sub(1);
        self.cursor.col = last_visible.min(self.active_alignment().alignment_length().saturating_sub(1));
    }

    /// Goes to a specific column (1-indexed for user).
    pub fn goto_column(&mut self, col: usize) {
        self.clear_pending();
        let aln_len = self.active_alignment().alignment_length();
        if col > 0 && col <= aln_len {
            self.cursor.col = col - 1; // 1-indexed for user
            self.ensure_cursor_visible();
            self.status_message = None; // Clear any previous error
        } else if col == 0 {
            self.cursor.col = 0;
            self.ensure_cursor_visible();
            self.status_message = None; // Clear any previous error
        } else {
            self.status_message = Some(format!("Invalid column: {}", col));
        }
    }

    /// Handles a g-command key.
    pub fn handle_g_command(&mut self, c: char) {
        self.pending_g = false;
        match c {
            '0' => self.goto_first_visible_column(),
            'm' => self.goto_middle_visible_column(),
            '$' => self.goto_last_visible_column(),
            _ => {} // Unknown g-command, ignore
        }
    }

    /// Cancels command mode and returns to normal mode.
    pub fn cancel_command(&mut self) {
        self.mode = AppMode::Normal;
    }

    /// Enters search mode.
    pub fn enter_search_mode(&mut self, backward: bool) {
        if backward {
            self.mode = AppMode::SearchBackward(String::new());
        } else {
            self.mode = AppMode::Search(String::new());
        }
    }

    /// Handles a character input in search mode.
    pub fn search_input(&mut self, c: char) {
        match &mut self.mode {
            AppMode::Search(ref mut pattern) | AppMode::SearchBackward(ref mut pattern) => {
                pattern.push(c);
            }
            _ => {}
        }
    }

    /// Handles backspace in search mode.
    pub fn search_backspace(&mut self) {
        match &mut self.mode {
            AppMode::Search(ref mut pattern) | AppMode::SearchBackward(ref mut pattern) => {
                pattern.pop();
                if pattern.is_empty() {
                    self.mode = AppMode::Normal;
                }
            }
            _ => {}
        }
    }

    /// Executes the current search.
    pub fn execute_search(&mut self) {
        let (pattern, backward) = match &self.mode {
            AppMode::Search(p) => (p.clone(), false),
            AppMode::SearchBackward(p) => (p.clone(), true),
            _ => return,
        };

        if pattern.is_empty() {
            self.mode = AppMode::Normal;
            return;
        }

        self.last_search = Some(pattern.clone());
        self.last_search_backward = backward;
        self.mode = AppMode::Normal;

        // Perform the search
        if backward {
            self.search_backward(&pattern);
        } else {
            self.search_forward(&pattern);
        }
    }

    /// Cancels search mode and returns to normal mode.
    pub fn cancel_search(&mut self) {
        self.mode = AppMode::Normal;
    }

    /// Finds the next match (n key).
    pub fn find_next(&mut self) {
        if let Some(pattern) = self.last_search.clone() {
            if self.last_search_backward {
                self.search_backward(&pattern);
            } else {
                self.search_forward(&pattern);
            }
        } else {
            self.status_message = Some("No previous search pattern".to_string());
        }
    }

    /// Finds the previous match (N key).
    pub fn find_previous(&mut self) {
        if let Some(pattern) = self.last_search.clone() {
            // Reverse the direction
            if self.last_search_backward {
                self.search_forward(&pattern);
            } else {
                self.search_backward(&pattern);
            }
        } else {
            self.status_message = Some("No previous search pattern".to_string());
        }
    }

    /// Searches forward from current position (right then down).
    /// Searches both sequence data and sequence names.
    fn search_forward(&mut self, pattern: &str) {
        let pattern_upper = pattern.to_uppercase();
        let start_row = self.cursor.row;
        let start_col = self.cursor.col;
        let alignment = self.active_alignment();
        let num_rows = alignment.sequence_count();
        let num_cols = alignment.alignment_length();

        if num_rows == 0 || num_cols == 0 {
            self.status_message = Some(format!("Pattern not found: {}", pattern));
            return;
        }

        // Search from current position to end of current row (sequence data only)
        let alignment = self.active_alignment();
        if let Some(seq) = alignment.get(start_row) {
            let search_start = start_col + 1; // Start after current position
            if search_start < seq.len() {
                let seq_str = seq.as_str();
                if let Some(pos) = seq_str[search_start..].to_uppercase().find(&pattern_upper) {
                    self.cursor.col = search_start + pos;
                    self.ensure_cursor_visible();
                    self.status_message = Some(format!("/{}", pattern));
                    return;
                }
            }
        }

        // Search remaining rows (name first, then data)
        for row_offset in 1..=num_rows {
            let row = (start_row + row_offset) % num_rows;
            let alignment = self.active_alignment();
            if let Some(seq) = alignment.get(row) {
                // First check sequence name
                if seq.id.to_uppercase().contains(&pattern_upper) {
                    self.cursor.row = row;
                    self.cursor.col = 0; // Position at start when matching name
                    self.ensure_cursor_visible();
                    let wrapped = row < start_row || (row == start_row);
                    if wrapped {
                        self.status_message = Some(format!("/{} (name, wrapped)", pattern));
                    } else {
                        self.status_message = Some(format!("/{} (name)", pattern));
                    }
                    return;
                }
                // Then check sequence data
                if let Some(pos) = seq.as_str().to_uppercase().find(&pattern_upper) {
                    self.cursor.row = row;
                    self.cursor.col = pos;
                    self.ensure_cursor_visible();
                    if row < start_row || (row == start_row && pos <= start_col) {
                        self.status_message = Some(format!("/{} (wrapped)", pattern));
                    } else {
                        self.status_message = Some(format!("/{}", pattern));
                    }
                    return;
                }
            }
        }

        self.status_message = Some(format!("Pattern not found: {}", pattern));
    }

    /// Searches backward from current position (left then up).
    /// Searches both sequence data and sequence names.
    fn search_backward(&mut self, pattern: &str) {
        let pattern_upper = pattern.to_uppercase();
        let start_row = self.cursor.row;
        let start_col = self.cursor.col;
        let alignment = self.active_alignment();
        let num_rows = alignment.sequence_count();
        let num_cols = alignment.alignment_length();

        if num_rows == 0 || num_cols == 0 {
            self.status_message = Some(format!("Pattern not found: {}", pattern));
            return;
        }

        // Search from current position backward in current row (data only)
        let alignment = self.active_alignment();
        if let Some(seq) = alignment.get(start_row) {
            if start_col > 0 {
                let search_area = &seq.as_str()[..start_col];
                if let Some(pos) = search_area.to_uppercase().rfind(&pattern_upper) {
                    self.cursor.col = pos;
                    self.ensure_cursor_visible();
                    self.status_message = Some(format!("?{}", pattern));
                    return;
                }
            }
            // Check current row's name if we're past position 0
            if start_col > 0 && seq.id.to_uppercase().contains(&pattern_upper) {
                self.cursor.col = 0;
                self.ensure_cursor_visible();
                self.status_message = Some(format!("?{} (name)", pattern));
                return;
            }
        }

        // Search previous rows (data first from end, then name)
        for row_offset in 1..=num_rows {
            let row = (start_row + num_rows - row_offset) % num_rows;
            let alignment = self.active_alignment();
            if let Some(seq) = alignment.get(row) {
                // First check sequence data (from end)
                if let Some(pos) = seq.as_str().to_uppercase().rfind(&pattern_upper) {
                    self.cursor.row = row;
                    self.cursor.col = pos;
                    self.ensure_cursor_visible();
                    if row > start_row || (row == start_row && pos >= start_col) {
                        self.status_message = Some(format!("?{} (wrapped)", pattern));
                    } else {
                        self.status_message = Some(format!("?{}", pattern));
                    }
                    return;
                }
                // Then check sequence name
                if seq.id.to_uppercase().contains(&pattern_upper) {
                    self.cursor.row = row;
                    self.cursor.col = 0;
                    self.ensure_cursor_visible();
                    let wrapped = row > start_row || (row == start_row);
                    if wrapped {
                        self.status_message = Some(format!("?{} (name, wrapped)", pattern));
                    } else {
                        self.status_message = Some(format!("?{} (name)", pattern));
                    }
                    return;
                }
            }
        }

        self.status_message = Some(format!("Pattern not found: {}", pattern));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequence_creation() {
        let seq = Sequence::new("seq1", "ACGT");
        assert_eq!(seq.id, "seq1");
        assert_eq!(seq.as_str(), "ACGT");
        assert_eq!(seq.len(), 4);
    }

    #[test]
    fn test_sequence_char_at() {
        let seq = Sequence::new("seq1", "ACGT");
        assert_eq!(seq.char_at(0), Some('A'));
        assert_eq!(seq.char_at(3), Some('T'));
        assert_eq!(seq.char_at(4), None);
    }

    #[test]
    fn test_alignment_valid() {
        let seqs = vec![
            Sequence::new("seq1", "ACGT"),
            Sequence::new("seq2", "TGCA"),
        ];
        let alignment = Alignment::new(seqs);
        assert!(alignment.is_valid_alignment);
        assert!(alignment.warning.is_none());
        assert_eq!(alignment.alignment_length(), 4);
    }

    #[test]
    fn test_alignment_invalid() {
        let seqs = vec![
            Sequence::new("seq1", "ACGT"),
            Sequence::new("seq2", "TG"),
        ];
        let alignment = Alignment::new(seqs);
        assert!(!alignment.is_valid_alignment);
        assert!(alignment.warning.is_some());
    }

    #[test]
    fn test_viewport_range() {
        let vp = Viewport::new(10, 20);
        assert_eq!(vp.row_range(), 0..10);
        assert_eq!(vp.col_range(), 0..20);
    }

    #[test]
    fn test_sequence_type_nucleotide() {
        let seqs = vec![
            Sequence::new("seq1", "ACGTACGT"),
            Sequence::new("seq2", "TGCA-TGC"),
        ];
        let alignment = Alignment::new(seqs);
        assert!(alignment.sequence_type.is_nucleotide());
        assert!(alignment.sequence_type.nt_ratio > 0.8);
    }

    #[test]
    fn test_sequence_type_amino_acid() {
        let seqs = vec![
            Sequence::new("seq1", "MKFLILLFNILCLFPVLAADNHGVGPQGAS"),
            Sequence::new("seq2", "MKWVTFISLLFLFSSAYSRGVFRRDAHKSE"),
        ];
        let alignment = Alignment::new(seqs);
        assert!(alignment.sequence_type.is_amino_acid());
        assert!(alignment.sequence_type.is_likely_not_nucleotide());
    }

    #[test]
    fn test_sequence_type_with_gaps() {
        // Amino acid sequence with gaps should still be detected
        let seqs = vec![
            Sequence::new("seq1", "MKF---LILLFNILCLFPVL"),
            Sequence::new("seq2", "MKW---VTFISLLFLFSSAY"),
        ];
        let alignment = Alignment::new(seqs);
        assert!(alignment.sequence_type.is_amino_acid());
    }

    #[test]
    fn test_cursor_movement() {
        let seqs = vec![
            Sequence::new("seq1", "ACGTACGT"),
            Sequence::new("seq2", "TGCATGCA"),
            Sequence::new("seq3", "AAAAAAAA"),
        ];
        let alignment = Alignment::new(seqs);
        let mut state = AppState::new(alignment, "test".to_string());
        state.update_viewport_size(3, 4);

        // Initial position
        assert_eq!(state.cursor.row, 0);
        assert_eq!(state.cursor.col, 0);

        // Move down
        state.move_down();
        assert_eq!(state.cursor.row, 1);

        // Move right
        state.move_right();
        assert_eq!(state.cursor.col, 1);

        // Move up
        state.move_up();
        assert_eq!(state.cursor.row, 0);

        // Move left
        state.move_left();
        assert_eq!(state.cursor.col, 0);

        // Boundary: can't go past 0
        state.move_up();
        assert_eq!(state.cursor.row, 0);
        state.move_left();
        assert_eq!(state.cursor.col, 0);
    }

    #[test]
    fn test_help_command() {
        let seqs = vec![Sequence::new("seq1", "ACGT")];
        let alignment = Alignment::new(seqs);
        let mut state = AppState::new(alignment, "test".to_string());

        // Initially help is not shown
        assert!(!state.show_help);

        // Enter command mode and type 'h'
        state.enter_command_mode();
        state.command_input('h');
        state.execute_command();

        // Should toggle help overlay on
        assert!(state.show_help);

        // Mode should return to normal
        assert_eq!(state.mode, AppMode::Normal);

        // Dismiss help
        state.dismiss_help();
        assert!(!state.show_help);
    }

    #[test]
    fn test_search_forward() {
        let seqs = vec![
            Sequence::new("seq1", "ACGTACGT"),
            Sequence::new("seq2", "TGCATGCA"),
            Sequence::new("seq3", "AAAAGAAA"),
        ];
        let alignment = Alignment::new(seqs);
        let mut state = AppState::new(alignment, "test".to_string());
        state.update_viewport_size(3, 8);

        // Search for "TGC" starting from (0,0)
        state.enter_search_mode(false);
        state.search_input('T');
        state.search_input('G');
        state.search_input('C');
        state.execute_search();

        // Should find in second row
        assert_eq!(state.cursor.row, 1);
        assert_eq!(state.cursor.col, 0);

        // Search for "GAA" - should find in third row
        state.enter_search_mode(false);
        state.search_input('G');
        state.search_input('A');
        state.search_input('A');
        state.execute_search();

        assert_eq!(state.cursor.row, 2);
        assert_eq!(state.cursor.col, 4);
    }

    #[test]
    fn test_search_backward() {
        let seqs = vec![
            Sequence::new("seq1", "ACGTACGT"),
            Sequence::new("seq2", "TGCATGCA"),
            Sequence::new("seq3", "AAAAGAAA"),
        ];
        let alignment = Alignment::new(seqs);
        let mut state = AppState::new(alignment, "test".to_string());
        state.update_viewport_size(3, 8);

        // Move to last row, last column
        state.cursor.row = 2;
        state.cursor.col = 7;

        // Search backward for "TGC"
        state.enter_search_mode(true);
        state.search_input('T');
        state.search_input('G');
        state.search_input('C');
        state.execute_search();

        // Should find last occurrence in second row (at position 4)
        assert_eq!(state.cursor.row, 1);
        assert_eq!(state.cursor.col, 4);
    }

    #[test]
    fn test_search_case_insensitive() {
        let seqs = vec![
            Sequence::new("seq1", "ACGTACGT"),
            Sequence::new("seq2", "tgcatgca"),
        ];
        let alignment = Alignment::new(seqs);
        let mut state = AppState::new(alignment, "test".to_string());
        state.update_viewport_size(2, 8);

        // Search for lowercase pattern in uppercase sequence
        state.enter_search_mode(false);
        state.search_input('c');
        state.search_input('g');
        state.search_input('t');
        state.execute_search();

        assert_eq!(state.cursor.row, 0);
        assert_eq!(state.cursor.col, 1);
    }

    #[test]
    fn test_find_next_and_previous() {
        let seqs = vec![
            Sequence::new("seq1", "ACGTACGT"),
            Sequence::new("seq2", "ACGTACGT"),
        ];
        let alignment = Alignment::new(seqs);
        let mut state = AppState::new(alignment, "test".to_string());
        state.update_viewport_size(2, 8);

        // Search for "CGT"
        state.enter_search_mode(false);
        state.search_input('C');
        state.search_input('G');
        state.search_input('T');
        state.execute_search();

        // First match at (0, 1)
        assert_eq!(state.cursor.row, 0);
        assert_eq!(state.cursor.col, 1);

        // Find next - should find at (0, 5)
        state.find_next();
        assert_eq!(state.cursor.row, 0);
        assert_eq!(state.cursor.col, 5);

        // Find next - should find in second row at (1, 1)
        state.find_next();
        assert_eq!(state.cursor.row, 1);
        assert_eq!(state.cursor.col, 1);

        // Find previous - should go back to (0, 5)
        state.find_previous();
        assert_eq!(state.cursor.row, 0);
        assert_eq!(state.cursor.col, 5);
    }

    #[test]
    fn test_search_not_found() {
        let seqs = vec![
            Sequence::new("seq1", "ACGTACGT"),
        ];
        let alignment = Alignment::new(seqs);
        let mut state = AppState::new(alignment, "test".to_string());
        state.update_viewport_size(1, 8);

        state.enter_search_mode(false);
        state.search_input('X');
        state.search_input('Y');
        state.search_input('Z');
        state.execute_search();

        // Cursor should not move
        assert_eq!(state.cursor.row, 0);
        assert_eq!(state.cursor.col, 0);
        assert!(state.status_message.as_ref().unwrap().contains("not found"));
    }

    #[test]
    fn test_search_by_sequence_name() {
        let seqs = vec![
            Sequence::new("Human_BRCA1", "ACGTACGT"),
            Sequence::new("Mouse_BRCA1", "TGCATGCA"),
            Sequence::new("Chicken_BRCA2", "AAAAGAAA"),
        ];
        let alignment = Alignment::new(seqs);
        let mut state = AppState::new(alignment, "test".to_string());
        state.update_viewport_size(3, 8);

        // Search for "Mouse" - should find by sequence name
        state.enter_search_mode(false);
        state.search_input('M');
        state.search_input('o');
        state.search_input('u');
        state.search_input('s');
        state.search_input('e');
        state.execute_search();

        assert_eq!(state.cursor.row, 1);
        assert_eq!(state.cursor.col, 0);
        assert!(state.status_message.as_ref().unwrap().contains("name"));

        // Search for "BRCA2" - should find Chicken sequence
        state.enter_search_mode(false);
        state.search_input('B');
        state.search_input('R');
        state.search_input('C');
        state.search_input('A');
        state.search_input('2');
        state.execute_search();

        assert_eq!(state.cursor.row, 2);
        assert_eq!(state.cursor.col, 0);
        assert!(state.status_message.as_ref().unwrap().contains("name"));
    }

    #[test]
    fn test_search_name_case_insensitive() {
        let seqs = vec![
            Sequence::new("Human_Gene", "ACGTACGT"),
            Sequence::new("mouse_gene", "TGCATGCA"),
        ];
        let alignment = Alignment::new(seqs);
        let mut state = AppState::new(alignment, "test".to_string());
        state.update_viewport_size(2, 8);

        // Search for lowercase "human" should match "Human_Gene"
        state.enter_search_mode(false);
        state.search_input('h');
        state.search_input('u');
        state.search_input('m');
        state.search_input('a');
        state.search_input('n');
        state.execute_search();

        // Should wrap around and find Human (row 0)
        assert_eq!(state.cursor.row, 0);
        assert!(state.status_message.as_ref().unwrap().contains("name"));
    }
}
