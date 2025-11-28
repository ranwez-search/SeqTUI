//! TUI rendering module.
//!
//! This module handles all visual rendering using ratatui:
//! - Layout with sticky sequence names on the left
//! - Colored nucleotide/amino acid display
//! - Status bar with position and mode info
//! - Hint bar with basic commands
//! - Help overlay popup
//!
//! The design supports future extensions like:
//! - Amino acid color schemes
//! - Codon highlighting
//! - Multiple panels (file browser, etc.)

use ratatui::{
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, Paragraph},
    Frame,
};

use crate::model::{AppMode, AppState};

/// Width reserved for sequence names (including border and padding).
const NAME_PANEL_WIDTH: u16 = 20;
/// Minimum width for the sequence panel.
const MIN_SEQ_PANEL_WIDTH: u16 = 10;
/// Height of the status bar.
const STATUS_BAR_HEIGHT: u16 = 1;
/// Height of the hint bar.
const HINT_BAR_HEIGHT: u16 = 1;

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
            'T' | 'U' => Color::Blue,
            // Gap - Light gray background
            '-' | '.' => Color::Rgb(180, 180, 180),
            // Unknown/ambiguous (N, etc.) - Medium gray
            'N' | 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' => Color::Rgb(140, 140, 140),
            // Stop/frameshift markers - White (must pop up!)
            '*' | '!' => Color::White,
            _ => Color::Rgb(140, 140, 140),
        }
    }
}

/// Amino acid color scheme based on Seaview's coloring.
/// Colors are chosen to match the familiar Seaview display.
pub struct AminoAcidColorScheme;

impl ColorScheme for AminoAcidColorScheme {
    fn get_color(&self, c: char) -> Color {
        match c.to_ascii_uppercase() {
            // Positively charged (basic) - Red
            'K' | 'R' => Color::Red,
            // Hydrophobic (nonpolar) - Blue
            'A' | 'F' | 'I' | 'L' | 'M' | 'V' | 'W' => Color::Blue,
            // Polar uncharged (green tones)
            'N' | 'Q' | 'S' | 'T' => Color::Green,
            // Aromatic Tyrosine and Histidine - Cyan
            'H' | 'Y' => Color::Cyan,
            // Cysteine - Pink/Light red
            'C' => Color::Rgb(255, 180, 180),
            // Negatively charged (acidic) - Magenta/Purple
            'D' | 'E' => Color::Magenta,
            // Proline - Yellow
            'P' => Color::Yellow,
            // Glycine - Light orange/salmon
            'G' => Color::Rgb(255, 200, 150),
            // Gap - Light gray background
            '-' | '.' => Color::Rgb(180, 180, 180),
            // Unknown/ambiguous (X, B, Z, etc.) - Medium gray (darker than gaps)
            'X' | 'B' | 'Z' | 'J' | 'O' | 'U' => Color::Rgb(140, 140, 140),
            // Stop codon and frameshift - White (must pop up!)
            '*' | '!' => Color::White,
            // Any other character - Medium gray
            _ => Color::Rgb(140, 140, 140),
        }
    }
}

/// Returns the appropriate color for a character based on sequence type.
fn get_color_for_sequence_type(c: char, seq_type: crate::model::SequenceType) -> Color {
    match seq_type {
        crate::model::SequenceType::Nucleotide => DnaColorScheme.get_color(c),
        crate::model::SequenceType::AminoAcid => AminoAcidColorScheme.get_color(c),
    }
}

/// Renders the complete UI.
pub fn render(frame: &mut Frame, state: &AppState) {
    let area = frame.area();

    // Main layout: content area + status bar + hint bar
    let main_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Min(3),
            Constraint::Length(STATUS_BAR_HEIGHT),
            Constraint::Length(HINT_BAR_HEIGHT),
        ])
        .split(area);

    let content_area = main_layout[0];
    let status_area = main_layout[1];
    let hint_area = main_layout[2];

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
    render_hint_bar(frame, hint_area);

    // Render help overlay if active
    if state.show_help {
        render_help_overlay(frame, area);
    }
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
    let seq_type = state.alignment.sequence_type;
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

                let bg_color = get_color_for_sequence_type(c, seq_type);
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

    // Show cursor position, sequence type, and visible range in title
    let type_str = match seq_type {
        crate::model::SequenceType::Nucleotide => "NT",
        crate::model::SequenceType::AminoAcid => "AA",
    };
    let title = format!(
        "Alignment ({}) [Site: {} | View: {}-{}/{}]",
        type_str,
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
        AppMode::Search(pattern) => ("SEARCH", format!("/{}", pattern)),
        AppMode::SearchBackward(pattern) => ("SEARCH", format!("?{}", pattern)),
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

/// Renders the hint bar at the very bottom with basic commands.
fn render_hint_bar(frame: &mut Frame, area: Rect) {
    let hints = vec![
        Span::styled(" :q ", Style::default().fg(Color::Black).bg(Color::DarkGray)),
        Span::styled(" quit ", Style::default().fg(Color::Gray)),
        Span::styled(" :h ", Style::default().fg(Color::Black).bg(Color::DarkGray)),
        Span::styled(" help ", Style::default().fg(Color::Gray)),
        Span::styled(" / ", Style::default().fg(Color::Black).bg(Color::DarkGray)),
        Span::styled(" search ", Style::default().fg(Color::Gray)),
        Span::styled(" h/j/k/l ", Style::default().fg(Color::Black).bg(Color::DarkGray)),
        Span::styled(" navigate ", Style::default().fg(Color::Gray)),
    ];

    let hint_line = Line::from(hints);
    let paragraph = Paragraph::new(hint_line);
    frame.render_widget(paragraph, area);
}

/// Renders a centered help overlay popup.
fn render_help_overlay(frame: &mut Frame, area: Rect) {
    // Calculate centered popup dimensions
    let popup_width = 60.min(area.width.saturating_sub(4));
    let popup_height = 22.min(area.height.saturating_sub(4));
    
    let popup_x = (area.width.saturating_sub(popup_width)) / 2;
    let popup_y = (area.height.saturating_sub(popup_height)) / 2;
    
    let popup_area = Rect::new(popup_x, popup_y, popup_width, popup_height);

    // Clear the area behind the popup
    frame.render_widget(Clear, popup_area);

    // Help content
    let help_lines = vec![
        Line::from(Span::styled("NAVIGATION", Style::default().add_modifier(Modifier::BOLD))),
        Line::from("  h / ←       Move left          l / →  Move right"),
        Line::from("  k / ↑       Move up            j / ↓  Move down"),
        Line::from("  0 / Home    First column       $ / End  Last column"),
        Line::from("  g0          First visible col  g$  Last visible col"),
        Line::from("  gm          Middle visible     gg  First column"),
        Line::from("  <num>|      Go to column (e.g., 50|)"),
        Line::from(""),
        Line::from(Span::styled("SEARCH", Style::default().add_modifier(Modifier::BOLD))),
        Line::from("  /pattern    Search forward (names & sequences)"),
        Line::from("  ?pattern    Search backward"),
        Line::from("  n           Next match         N   Previous match"),
        Line::from(""),
        Line::from(Span::styled("COMMANDS", Style::default().add_modifier(Modifier::BOLD))),
        Line::from("  :q          Quit               :h  Toggle help"),
        Line::from("  :<number>   Go to column"),
        Line::from(""),
        Line::from(Span::styled("Press any key to close", Style::default().fg(Color::DarkGray))),
    ];

    let block = Block::default()
        .borders(Borders::ALL)
        .title(" Help ")
        .title_style(Style::default().add_modifier(Modifier::BOLD))
        .style(Style::default().bg(Color::Black));

    let paragraph = Paragraph::new(help_lines).block(block);
    frame.render_widget(paragraph, popup_area);
}

/// Calculates the visible dimensions for the sequence panel.
pub fn calculate_visible_dimensions(terminal_width: u16, terminal_height: u16) -> (usize, usize) {
    // Account for name panel, sequence panel borders, status bar, and hint bar
    // Sequence panel width = terminal_width - NAME_PANEL_WIDTH
    // Visible cols = sequence panel width - 2 (for left/right borders)
    let visible_cols = (terminal_width.saturating_sub(NAME_PANEL_WIDTH + 2)) as usize;
    let visible_rows = (terminal_height.saturating_sub(STATUS_BAR_HEIGHT + HINT_BAR_HEIGHT + 2)) as usize;
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
        assert_eq!(scheme.get_color('-'), Color::Rgb(180, 180, 180)); // Light gray for gaps
        assert_eq!(scheme.get_color('N'), Color::Rgb(140, 140, 140)); // Medium gray for unknown
        assert_eq!(scheme.get_color('*'), Color::White); // Stop codon pops up
    }

    #[test]
    fn test_amino_acid_colors() {
        let scheme = AminoAcidColorScheme;
        assert_eq!(scheme.get_color('K'), Color::Red);     // Basic
        assert_eq!(scheme.get_color('L'), Color::Blue);    // Hydrophobic
        assert_eq!(scheme.get_color('S'), Color::Green);   // Polar
        assert_eq!(scheme.get_color('P'), Color::Yellow);  // Proline
        assert_eq!(scheme.get_color('-'), Color::Rgb(180, 180, 180)); // Light gray for gaps
        assert_eq!(scheme.get_color('X'), Color::Rgb(140, 140, 140)); // Medium gray for unknown
        assert_eq!(scheme.get_color('*'), Color::White);   // Stop codon pops up
    }

    #[test]
    fn test_visible_dimensions() {
        let (rows, cols) = calculate_visible_dimensions(100, 50);
        // 100 - 20 (name panel) - 2 (sequence panel borders) = 78 cols
        // 50 - 1 (status) - 1 (hint) - 2 (borders) = 46 rows
        assert_eq!(cols, 78);
        assert_eq!(rows, 46);
    }
}
