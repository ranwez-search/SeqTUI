//! Data model for the alignment viewer.
//!
//! This module contains all data structures for representing:
//! - Sequences and alignments
//! - Viewport state
//! - Application state
//!
//! The design allows for future extensions like filtering, codon views,
//! and translation between nucleotides and amino acids.

use std::ops::Range;

/// Represents a single sequence with its identifier and data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    /// The sequence identifier (from FASTA header, without '>')
    pub id: String,
    /// The sequence data (nucleotides or amino acids)
    pub data: String,
}

impl Sequence {
    /// Creates a new sequence.
    pub fn new(id: impl Into<String>, data: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            data: data.into(),
        }
    }

    /// Returns the length of the sequence.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Returns true if the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Gets a character at a specific position.
    /// 
    /// This uses byte indexing for O(1) performance, which is safe for
    /// DNA/amino acid sequences that are always ASCII.
    pub fn char_at(&self, pos: usize) -> Option<char> {
        // Use byte indexing for O(1) access - safe for ASCII sequences
        self.data.as_bytes().get(pos).map(|&b| b as char)
    }

    /// Gets a slice of the sequence data.
    pub fn slice(&self, range: Range<usize>) -> &str {
        let start = range.start.min(self.data.len());
        let end = range.end.min(self.data.len());
        &self.data[start..end]
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
}

impl Alignment {
    /// Creates a new alignment from a vector of sequences.
    pub fn new(sequences: Vec<Sequence>) -> Self {
        let (is_valid, alignment_length, warning) = Self::validate_alignment(&sequences);
        Self {
            sequences,
            alignment_length,
            is_valid_alignment: is_valid,
            warning,
        }
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
}

/// The complete application state.
#[derive(Debug)]
pub struct AppState {
    /// The loaded alignment
    pub alignment: Alignment,
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
    /// Pending 'g' key for g-commands
    pub pending_g: bool,
    /// Number buffer for <number>| command
    pub number_buffer: String,
}

impl AppState {
    /// Creates a new application state with the given alignment.
    pub fn new(alignment: Alignment) -> Self {
        let warning = alignment.warning.clone();
        Self {
            alignment,
            viewport: Viewport::new(0, 0),
            cursor: Cursor::new(),
            mode: AppMode::Normal,
            should_quit: false,
            status_message: warning,
            last_search: None,
            last_search_backward: false,
            show_help: false,
            pending_g: false,
            number_buffer: String::new(),
        }
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
        if self.cursor.row + 1 < self.alignment.sequence_count() {
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
        if self.cursor.col + 1 < self.alignment.alignment_length() {
            self.cursor.col += 1;
            self.ensure_cursor_visible();
        }
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
        let max_row = self.alignment.sequence_count().saturating_sub(1);
        let max_col = self.alignment.alignment_length().saturating_sub(1);

        // Ensure we don't scroll past the end
        if self.viewport.first_row + self.viewport.visible_rows > self.alignment.sequence_count() {
            self.viewport.first_row = self.alignment.sequence_count().saturating_sub(self.viewport.visible_rows);
        }

        if self.viewport.first_col + self.viewport.visible_cols > self.alignment.alignment_length() {
            self.viewport.first_col = self.alignment.alignment_length().saturating_sub(self.viewport.visible_cols);
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
    pub fn execute_command(&mut self) {
        if let AppMode::Command(ref cmd) = self.mode.clone() {
            match cmd.as_str() {
                "q" | "quit" => self.should_quit = true,
                "h" | "help" => self.show_help(),
                _ => {
                    // Future: handle :number for column navigation
                    if let Ok(col) = cmd.parse::<usize>() {
                        if col > 0 && col <= self.alignment.alignment_length() {
                            self.cursor.col = col - 1; // 1-indexed for user
                            self.ensure_cursor_visible();
                        } else {
                            self.status_message = Some(format!("Invalid column: {}", col));
                        }
                    } else {
                        self.status_message = Some(format!("Unknown command: {}", cmd));
                    }
                }
            }
        }
        self.mode = AppMode::Normal;
    }

    /// Shows help information.
    fn show_help(&mut self) {
        self.show_help = !self.show_help;
    }

    /// Toggles help off (for any key press while help is shown).
    pub fn dismiss_help(&mut self) {
        self.show_help = false;
    }

    /// Clears any pending key state.
    fn clear_pending(&mut self) {
        self.pending_g = false;
        self.number_buffer.clear();
    }

    /// Sets the pending 'g' state for g-commands.
    pub fn set_pending_g(&mut self) {
        self.pending_g = true;
        self.number_buffer.clear();
    }

    /// Accumulates a digit for the number prefix.
    pub fn accumulate_digit(&mut self, c: char) {
        self.pending_g = false;
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
        if self.alignment.alignment_length() > 0 {
            self.cursor.col = self.alignment.alignment_length() - 1;
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
        self.cursor.col = middle.min(self.alignment.alignment_length().saturating_sub(1));
    }

    /// Goes to the last visible column (g$).
    pub fn goto_last_visible_column(&mut self) {
        self.clear_pending();
        let last_visible = self.viewport.first_col + self.viewport.visible_cols.saturating_sub(1);
        self.cursor.col = last_visible.min(self.alignment.alignment_length().saturating_sub(1));
    }

    /// Goes to a specific column (1-indexed for user).
    pub fn goto_column(&mut self, col: usize) {
        self.clear_pending();
        if col > 0 && col <= self.alignment.alignment_length() {
            self.cursor.col = col - 1; // 1-indexed for user
            self.ensure_cursor_visible();
        } else if col == 0 {
            self.cursor.col = 0;
            self.ensure_cursor_visible();
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
        let num_rows = self.alignment.sequence_count();
        let num_cols = self.alignment.alignment_length();

        if num_rows == 0 || num_cols == 0 {
            self.status_message = Some(format!("Pattern not found: {}", pattern));
            return;
        }

        // Search from current position to end of current row (sequence data only)
        if let Some(seq) = self.alignment.get(start_row) {
            let search_start = start_col + 1; // Start after current position
            if search_start < seq.data.len() {
                if let Some(pos) = seq.data[search_start..].to_uppercase().find(&pattern_upper) {
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
            if let Some(seq) = self.alignment.get(row) {
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
                if let Some(pos) = seq.data.to_uppercase().find(&pattern_upper) {
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
        let num_rows = self.alignment.sequence_count();
        let num_cols = self.alignment.alignment_length();

        if num_rows == 0 || num_cols == 0 {
            self.status_message = Some(format!("Pattern not found: {}", pattern));
            return;
        }

        // Search from current position backward in current row (data only)
        if let Some(seq) = self.alignment.get(start_row) {
            if start_col > 0 {
                let search_area = &seq.data[..start_col];
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
            if let Some(seq) = self.alignment.get(row) {
                // First check sequence data (from end)
                if let Some(pos) = seq.data.to_uppercase().rfind(&pattern_upper) {
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
        assert_eq!(seq.data, "ACGT");
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
    fn test_cursor_movement() {
        let seqs = vec![
            Sequence::new("seq1", "ACGTACGT"),
            Sequence::new("seq2", "TGCATGCA"),
            Sequence::new("seq3", "AAAAAAAA"),
        ];
        let alignment = Alignment::new(seqs);
        let mut state = AppState::new(alignment);
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
        let mut state = AppState::new(alignment);

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
        let mut state = AppState::new(alignment);
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
        let mut state = AppState::new(alignment);
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
        let mut state = AppState::new(alignment);
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
        let mut state = AppState::new(alignment);
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
        let mut state = AppState::new(alignment);
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
        let mut state = AppState::new(alignment);
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
        let mut state = AppState::new(alignment);
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
