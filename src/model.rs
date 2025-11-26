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
    pub fn char_at(&self, pos: usize) -> Option<char> {
        self.data.chars().nth(pos)
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

    /// Cancels command mode and returns to normal mode.
    pub fn cancel_command(&mut self) {
        self.mode = AppMode::Normal;
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
}
