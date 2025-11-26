//! TUI rendering module.
//!
//! This module handles all visual rendering using ratatui:
//! - Layout with sticky sequence names on the left
//! - Colored nucleotide/amino acid display
//! - Status bar with position and mode info
//! - Command line display
//!
//! The design supports future extensions like:
//! - Amino acid color schemes
//! - Codon highlighting
//! - Multiple panels (file browser, etc.)

use ratatui::{
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Paragraph},
    Frame,
};

use crate::model::{AppMode, AppState};

/// Width reserved for sequence names (including border and padding).
const NAME_PANEL_WIDTH: u16 = 20;
/// Minimum width for the sequence panel.
const MIN_SEQ_PANEL_WIDTH: u16 = 10;
/// Height of the status bar.
const STATUS_BAR_HEIGHT: u16 = 1;

/// Color scheme for nucleotides.
///
/// This trait allows for different color schemes to be implemented
/// (e.g., for amino acids in the future).
pub trait ColorScheme {
    fn get_color(&self, c: char) -> Color;
}

/// DNA nucleotide color scheme.
pub struct DnaColorScheme;

impl ColorScheme for DnaColorScheme {
    fn get_color(&self, c: char) -> Color {
        match c.to_ascii_uppercase() {
            'A' => Color::Red,
            'C' => Color::Green,
            'G' => Color::Yellow,
            'T' => Color::Blue,
            _ => Color::DarkGray,
        }
    }
}

/// Amino acid color scheme (placeholder for future implementation).
#[allow(dead_code)]
pub struct AminoAcidColorScheme;

impl ColorScheme for AminoAcidColorScheme {
    fn get_color(&self, c: char) -> Color {
        // Placeholder - implement proper amino acid colors later
        // Could use physicochemical properties for grouping
        match c.to_ascii_uppercase() {
            // Hydrophobic
            'A' | 'V' | 'I' | 'L' | 'M' | 'F' | 'W' | 'P' => Color::Yellow,
            // Polar
            'S' | 'T' | 'N' | 'Q' | 'C' | 'G' | 'Y' => Color::Green,
            // Charged positive
            'K' | 'R' | 'H' => Color::Blue,
            // Charged negative
            'D' | 'E' => Color::Red,
            // Gap or unknown
            '-' | 'X' | '*' => Color::DarkGray,
            _ => Color::Gray,
        }
    }
}

/// Renders the complete UI.
pub fn render(frame: &mut Frame, state: &AppState) {
    let area = frame.area();

    // Main layout: content area + status bar
    let main_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Min(3),
            Constraint::Length(STATUS_BAR_HEIGHT),
        ])
        .split(area);

    let content_area = main_layout[0];
    let status_area = main_layout[1];

    // Split content area: names panel (left) + sequence panel (right)
    let content_layout = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Length(NAME_PANEL_WIDTH),
            Constraint::Min(MIN_SEQ_PANEL_WIDTH),
        ])
        .split(content_area);

    let names_area = content_layout[0];
    let sequences_area = content_layout[1];

    // Calculate visible dimensions (accounting for borders)
    let visible_rows = (sequences_area.height.saturating_sub(2)) as usize; // -2 for borders
    let visible_cols = (sequences_area.width.saturating_sub(2)) as usize; // -2 for borders

    // Update state viewport if dimensions changed
    // (This is a bit of a hack - ideally state update would be in controller)
    // We'll handle this properly in the main loop

    // Render each panel
    render_names_panel(frame, state, names_area, visible_rows);
    render_sequences_panel(frame, state, sequences_area, visible_rows, visible_cols);
    render_status_bar(frame, state, status_area);
}

/// Renders the sequence names panel (sticky, always visible).
fn render_names_panel(frame: &mut Frame, state: &AppState, area: Rect, visible_rows: usize) {
    let mut lines: Vec<Line> = Vec::new();

    let start_row = state.viewport.first_row;
    let end_row = (start_row + visible_rows).min(state.alignment.sequence_count());

    for row_idx in start_row..end_row {
        if let Some(seq) = state.alignment.get(row_idx) {
            let is_current = row_idx == state.cursor.row;

            // Truncate name if too long
            let max_name_len = (NAME_PANEL_WIDTH.saturating_sub(3)) as usize;
            let name = if seq.id.len() > max_name_len {
                format!("{}…", &seq.id[..max_name_len - 1])
            } else {
                seq.id.clone()
            };

            let style = if is_current {
                Style::default()
                    .fg(Color::Black)
                    .bg(Color::White)
                    .add_modifier(Modifier::BOLD)
            } else {
                Style::default().fg(Color::White)
            };

            lines.push(Line::from(Span::styled(name, style)));
        }
    }

    let block = Block::default()
        .borders(Borders::ALL)
        .title("Sequences");

    let paragraph = Paragraph::new(lines).block(block);
    frame.render_widget(paragraph, area);
}

/// Renders the sequences panel with colored nucleotides/amino acids.
fn render_sequences_panel(
    frame: &mut Frame,
    state: &AppState,
    area: Rect,
    visible_rows: usize,
    visible_cols: usize,
) {
    let color_scheme = DnaColorScheme;
    let mut lines: Vec<Line> = Vec::new();

    let start_row = state.viewport.first_row;
    let end_row = (start_row + visible_rows).min(state.alignment.sequence_count());
    let start_col = state.viewport.first_col;
    let end_col = (start_col + visible_cols).min(state.alignment.alignment_length());

    for row_idx in start_row..end_row {
        if let Some(seq) = state.alignment.get(row_idx) {
            let is_current_row = row_idx == state.cursor.row;
            let mut spans: Vec<Span> = Vec::new();

            for col_idx in start_col..end_col {
                let c = seq.char_at(col_idx).unwrap_or(' ');
                let is_cursor = is_current_row && col_idx == state.cursor.col;

                let bg_color = color_scheme.get_color(c);
                let fg_color = Color::Black;

                let style = if is_cursor {
                    // Invert colors for cursor position
                    Style::default()
                        .fg(bg_color)
                        .bg(Color::White)
                        .add_modifier(Modifier::BOLD)
                } else {
                    Style::default().fg(fg_color).bg(bg_color)
                };

                spans.push(Span::styled(c.to_string(), style));
            }

            lines.push(Line::from(spans));
        }
    }

    // Show cursor position and visible range in title
    let title = format!(
        "Alignment [Site: {} | View: {}-{}/{}]",
        state.cursor.col + 1,
        start_col + 1,
        end_col,
        state.alignment.alignment_length()
    );

    let block = Block::default().borders(Borders::ALL).title(title);

    let paragraph = Paragraph::new(lines).block(block);
    frame.render_widget(paragraph, area);
}

/// Renders the status bar at the bottom.
fn render_status_bar(frame: &mut Frame, state: &AppState, area: Rect) {
    let (mode_str, command_str) = match &state.mode {
        AppMode::Normal => ("NORMAL", String::new()),
        AppMode::Command(cmd) => ("COMMAND", format!(":{}", cmd)),
    };

    let position_info = format!(
        "Seq {}/{} | Col {}/{} ",
        state.cursor.row + 1,
        state.alignment.sequence_count(),
        state.cursor.col + 1,
        state.alignment.alignment_length()
    );

    // Show warning or status message if present
    let message = state.status_message.as_deref().unwrap_or("");

    let left_content = if command_str.is_empty() {
        format!(" {} | {} ", mode_str, message)
    } else {
        format!(" {} | {} ", mode_str, command_str)
    };

    let left_len = left_content.len();
    let status_line = Line::from(vec![
        Span::styled(
            left_content,
            Style::default().fg(Color::Black).bg(Color::Cyan),
        ),
        Span::styled(
            " ".repeat((area.width as usize).saturating_sub(left_len + position_info.len())),
            Style::default().bg(Color::Cyan),
        ),
        Span::styled(
            position_info,
            Style::default()
                .fg(Color::Black)
                .bg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
    ]);

    let paragraph = Paragraph::new(status_line);
    frame.render_widget(paragraph, area);
}

/// Calculates the visible dimensions for the sequence panel.
pub fn calculate_visible_dimensions(terminal_width: u16, terminal_height: u16) -> (usize, usize) {
    // Account for borders and status bar
    let visible_cols = (terminal_width.saturating_sub(NAME_PANEL_WIDTH + 4)) as usize;
    let visible_rows = (terminal_height.saturating_sub(STATUS_BAR_HEIGHT + 2)) as usize;
    (visible_rows, visible_cols)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dna_colors() {
        let scheme = DnaColorScheme;
        assert_eq!(scheme.get_color('A'), Color::Red);
        assert_eq!(scheme.get_color('a'), Color::Red); // Case insensitive
        assert_eq!(scheme.get_color('C'), Color::Green);
        assert_eq!(scheme.get_color('G'), Color::Yellow);
        assert_eq!(scheme.get_color('T'), Color::Blue);
        assert_eq!(scheme.get_color('-'), Color::DarkGray);
        assert_eq!(scheme.get_color('N'), Color::DarkGray);
    }

    #[test]
    fn test_visible_dimensions() {
        let (rows, cols) = calculate_visible_dimensions(100, 50);
        // 100 - 20 (name panel) - 4 (borders) = 76 cols
        // 50 - 1 (status) - 2 (borders) = 47 rows
        assert_eq!(cols, 76);
        assert_eq!(rows, 47);
    }
}
