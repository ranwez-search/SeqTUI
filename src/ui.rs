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

mod glyphs;

use ratatui::{
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, Paragraph},
    Frame,
};

use crate::model::{AppMode, AppState, ViewMode};
use glyphs::Glyphs;

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
    if seq_type.is_nucleotide() {
        DnaColorScheme.get_color(c)
    } else {
        AminoAcidColorScheme.get_color(c)
    }
}

/// Renders the complete UI.
pub fn render(frame: &mut Frame, state: &AppState) {
    let glyphs = glyphs::select(state.fancy_ui);
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
    render_status_bar(frame, state, &glyphs, status_area);
    render_hint_bar(frame, &glyphs, hint_area);

    // Render help overlay if active
    if state.show_help {
        render_help_overlay(frame, state, &glyphs, area);
    }

    // Render translation settings overlay if active
    if state.mode == AppMode::TranslationSettings {
        render_translation_settings_overlay(frame, state, &glyphs, area);
    }

    // Render loading overlay if active
    if state.loading_state.is_loading() {
        render_loading_overlay(frame, state, area);
    }

    // Render error popup if active
    if let Some(error_msg) = &state.error_popup {
        render_error_popup(frame, error_msg, area);
    }

    // Render file browser if active
    if let Some(browser) = &state.file_browser {
        render_file_browser(frame, browser, &glyphs, area);
}
}

/// Renders the sequence names panel (sticky, always visible).
fn render_names_panel(frame: &mut Frame, state: &AppState, area: Rect, visible_rows: usize) {
    let mut lines: Vec<Line> = Vec::new();
    let alignment = state.active_alignment();

    let start_row = state.viewport.first_row;
    let end_row = (start_row + visible_rows).min(alignment.sequence_count());

    for row_idx in start_row..end_row {
        if let Some(seq) = alignment.get(row_idx) {
            let is_current = row_idx == state.cursor.row;

            // Truncate name if too long
            let max_name_len = (NAME_PANEL_WIDTH.saturating_sub(3)) as usize;
            let name = if seq.id.len() > max_name_len {
                let truncate_len = max_name_len.saturating_sub(3);
                format!("{}...", &seq.id[..truncate_len])
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
    let alignment = state.active_alignment();
    let seq_type = alignment.sequence_type;
    let mut lines: Vec<Line> = Vec::new();

    let start_row = state.viewport.first_row;
    let end_row = (start_row + visible_rows).min(alignment.sequence_count());
    let start_col = state.viewport.first_col;
    let end_col = (start_col + visible_cols).min(alignment.alignment_length());

    for row_idx in start_row..end_row {
        if let Some(seq) = alignment.get(row_idx) {
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

    // Show file name, view mode info, and visible range in title
    let view_info = match state.view_mode {
        ViewMode::Nucleotide => "NT".to_string(),
        ViewMode::AminoAcid => format!(
            "AA, code {}, frame +{}",
            state.translation_settings.genetic_code_id,
            state.translation_settings.frame + 1
        ),
    };
    let title = format!(
        "{} ({}) [Site: {} | View: {}-{}/{}]",
        state.file_name,
        view_info,
        state.cursor.col + 1,
        start_col + 1,
        end_col,
        alignment.alignment_length()
    );

    let block = Block::default().borders(Borders::ALL).title(title);

    let paragraph = Paragraph::new(lines).block(block);
    frame.render_widget(paragraph, area);
}

/// Renders the status bar at the bottom.
fn render_status_bar(frame: &mut Frame, state: &AppState, glyphs: &Glyphs, area: Rect) {
    let alignment = state.active_alignment();
    
    let position_info = format!(
        "Seq {}/{} | Col {}/{} ",
        state.cursor.row + 1,
        alignment.sequence_count(),
        state.cursor.col + 1,
        alignment.alignment_length()
    );

    // Build status line based on mode
    let status_line = match &state.mode {
        AppMode::Normal => {
            // Show status message or just NORMAL
            let message = state.status_message.as_deref().unwrap_or("");
            let left_content = format!(" NORMAL | {} ", message);
            let left_len = left_content.len();
            
            Line::from(vec![
                Span::styled(
                    left_content,
                    Style::default().fg(Color::Black).bg(Color::Cyan),
                ),
                Span::styled(
                    " ".repeat((area.width as usize).saturating_sub(left_len + position_info.len())),
                    Style::default().bg(Color::Cyan),
                ),
                Span::styled(
                    position_info.clone(),
                    Style::default().fg(Color::Black).bg(Color::Cyan).add_modifier(Modifier::BOLD),
                ),
            ])
        }
        AppMode::Command(cmd) => {
            // Mode label with cyan background, input area with white background + cursor
            let mode_label = " COMMAND ";
            let input_with_cursor = format!(":{}{}", cmd, glyphs.cursor);
            let input_len = input_with_cursor.len() + 1; // +1 for space after
            let mode_len = mode_label.len();
            
            Line::from(vec![
                Span::styled(
                    mode_label,
                    Style::default().fg(Color::Black).bg(Color::Cyan).add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    format!(" {} ", input_with_cursor),
                    Style::default().fg(Color::Black).bg(Color::White),
                ),
                Span::styled(
                    " ".repeat((area.width as usize).saturating_sub(mode_len + input_len + 1 + position_info.len())),
                    Style::default().bg(Color::Cyan),
                ),
                Span::styled(
                    position_info.clone(),
                    Style::default().fg(Color::Black).bg(Color::Cyan).add_modifier(Modifier::BOLD),
                ),
            ])
        }
        AppMode::Search(pattern) => {
            // Search forward mode
            let mode_label = " SEARCH ";
            let input_with_cursor = format!("/{}{}", pattern, glyphs.cursor);
            let input_len = input_with_cursor.len() + 1;
            let mode_len = mode_label.len();
            
            Line::from(vec![
                Span::styled(
                    mode_label,
                    Style::default().fg(Color::Black).bg(Color::Cyan).add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    format!(" {} ", input_with_cursor),
                    Style::default().fg(Color::Black).bg(Color::White),
                ),
                Span::styled(
                    " ".repeat((area.width as usize).saturating_sub(mode_len + input_len + 1 + position_info.len())),
                    Style::default().bg(Color::Cyan),
                ),
                Span::styled(
                    position_info.clone(),
                    Style::default().fg(Color::Black).bg(Color::Cyan).add_modifier(Modifier::BOLD),
                ),
            ])
        }
        AppMode::SearchBackward(pattern) => {
            // Search backward mode
            let mode_label = " SEARCH ";
            let input_with_cursor = format!("?{}{}", pattern, glyphs.cursor);
            let input_len = input_with_cursor.len() + 1;
            let mode_len = mode_label.len();
            
            Line::from(vec![
                Span::styled(
                    mode_label,
                    Style::default().fg(Color::Black).bg(Color::Cyan).add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    format!(" {} ", input_with_cursor),
                    Style::default().fg(Color::Black).bg(Color::White),
                ),
                Span::styled(
                    " ".repeat((area.width as usize).saturating_sub(mode_len + input_len + 1 + position_info.len())),
                    Style::default().bg(Color::Cyan),
                ),
                Span::styled(
                    position_info.clone(),
                    Style::default().fg(Color::Black).bg(Color::Cyan).add_modifier(Modifier::BOLD),
                ),
            ])
        }
        AppMode::TranslationSettings => {
            let left_content = " TRANSLATE | Use j/k for code, h/l for frame ";
            let left_len = left_content.len();
            
            Line::from(vec![
                Span::styled(
                    left_content,
                    Style::default().fg(Color::Black).bg(Color::Cyan),
                ),
                Span::styled(
                    " ".repeat((area.width as usize).saturating_sub(left_len + position_info.len())),
                    Style::default().bg(Color::Cyan),
                ),
                Span::styled(
                    position_info.clone(),
                    Style::default().fg(Color::Black).bg(Color::Cyan).add_modifier(Modifier::BOLD),
                ),
            ])
        }
    };

    let paragraph = Paragraph::new(status_line);
    frame.render_widget(paragraph, area);
}

/// Renders the hint bar at the very bottom with basic commands.
fn render_hint_bar(frame: &mut Frame, glyphs: &Glyphs, area: Rect) {
    let arrows = format!(
        " {}{}{}{} ",
        glyphs.arrow_left, glyphs.arrow_up, glyphs.arrow_down, glyphs.arrow_right
    );
    let hints = vec![
        Span::styled(" :q ", Style::default().fg(Color::Black).bg(Color::DarkGray)),
        Span::styled(" quit ", Style::default().fg(Color::Gray)),
        Span::styled(" :h ", Style::default().fg(Color::Black).bg(Color::DarkGray)),
        Span::styled(" help ", Style::default().fg(Color::Gray)),
        Span::styled(" / ", Style::default().fg(Color::Black).bg(Color::DarkGray)),
        Span::styled(" search ", Style::default().fg(Color::Gray)),
        Span::styled(arrows, Style::default().fg(Color::Black).bg(Color::DarkGray)),
        Span::styled(" navigate ", Style::default().fg(Color::Gray)),
    ];

    let hint_line = Line::from(hints);
    let paragraph = Paragraph::new(hint_line);
    frame.render_widget(paragraph, area);
}

/// Renders a centered help overlay popup with tabs.
fn render_help_overlay(frame: &mut Frame, state: &AppState, glyphs: &Glyphs, area: Rect) {
    use crate::model::HelpTab;
    
    // Calculate centered popup dimensions
    let popup_width = 58.min(area.width.saturating_sub(4));
    let popup_height = 18.min(area.height.saturating_sub(4));
    
    let popup_x = (area.width.saturating_sub(popup_width)) / 2;
    let popup_y = (area.height.saturating_sub(popup_height)) / 2;
    
    let popup_area = Rect::new(popup_x, popup_y, popup_width, popup_height);

    // Clear the area behind the popup
    frame.render_widget(Clear, popup_area);

    // Build tab bar
    let tab_spans: Vec<Span> = HelpTab::all()
        .iter()
        .map(|&tab| {
            let name = format!(" {} ", tab.name());
            if tab == state.help_tab {
                Span::styled(
                    name,
                    Style::default()
                        .fg(Color::Black)
                        .bg(Color::Cyan)
                        .add_modifier(Modifier::BOLD),
                )
            } else {
                Span::styled(name, Style::default().fg(Color::Gray))
            }
        })
        .collect();

    // Help content based on selected tab
    let mut help_lines: Vec<Line> = vec![
        Line::from(tab_spans),
        Line::from(""),
    ];

    match state.help_tab {
        HelpTab::Basics => {
            let tab_hint = format!(
                "Use {}/{} or h/l to switch tabs",
                glyphs.arrow_left, glyphs.arrow_right
            );
            help_lines.extend(vec![
                Line::from(Span::styled("QUICK COMMANDS", Style::default().add_modifier(Modifier::BOLD))),
                Line::from(""),
                Line::from("  :q             Quit the application"),
                Line::from("  :h             Toggle this help"),
                Line::from("  :e             Open file browser"),
                Line::from("  :<number>      Jump to sequence/row number"),
                Line::from("  :w file.fa     Save to FASTA (single-line seqs)"),
                Line::from(""),
                Line::from(Span::styled("CLI MODE", Style::default().add_modifier(Modifier::BOLD))),
                Line::from(""),
                Line::from("  Run 'seqtui -h' for CLI options (convert,"),
                Line::from("  concatenate, translate, VCF export)."),
                Line::from("  https://github.com/ranwez-search/SeqTUI"),
                Line::from(""),
                Line::from(Span::styled(tab_hint, Style::default().fg(Color::DarkGray))),
            ]);
        }
        HelpTab::Navigation => {
            let arrows = format!(
                "{}{}{}{}",
                glyphs.arrow_left, glyphs.arrow_up, glyphs.arrow_down, glyphs.arrow_right
            );
            let lr = format!("{}{}", glyphs.arrow_left, glyphs.arrow_right);
            let ud = format!("{}{}", glyphs.arrow_up, glyphs.arrow_down);
            let tab_hint = format!(
                "Use {}/{} or h/l to switch tabs",
                glyphs.arrow_left, glyphs.arrow_right
            );
            help_lines.extend(vec![
                Line::from(Span::styled("ARROW KEY NAVIGATION", Style::default().add_modifier(Modifier::BOLD))),
                Line::from(""),
                Line::from(format!("  {}           Move one position", arrows)),
                Line::from(format!("  Shift + {}     Half page left/right", lr)),
                Line::from(format!("  Shift + {}     Full page up/down", ud)),
                Line::from(""),
                Line::from("  Home           Go to first column"),
                Line::from("  End            Go to last column"),
                Line::from("  PgUp / PgDn    Page up/down"),
                Line::from(""),
                Line::from(Span::styled(tab_hint, Style::default().fg(Color::DarkGray))),
            ]);
        }
        HelpTab::VimNav => {
            let tab_hint = format!(
                "Use {}/{} or h/l to switch tabs",
                glyphs.arrow_left, glyphs.arrow_right
            );
            help_lines.extend(vec![
                Line::from(Span::styled("VIM-STYLE NAVIGATION", Style::default().add_modifier(Modifier::BOLD))),
                Line::from(""),
                Line::from("  h / j / k / l  Move left/down/up/right"),
                Line::from("  w / b / e      Next/prev/end of word"),
                Line::from("  Ctrl+U / D     Half page up/down"),
                Line::from("  zH / zL        Half page left/right"),
                Line::from(""),
                Line::from("  0 / $          First / last column"),
                Line::from("  g0 / gm / g$   First/middle/last visible column"),
                Line::from("  <num>|         Go to column (e.g., 50|)"),
                Line::from(""),
                Line::from(Span::styled(tab_hint, Style::default().fg(Color::DarkGray))),
            ]);
        }
        HelpTab::Search => {
            let tab_hint = format!(
                "Use {}/{} or h/l to switch tabs",
                glyphs.arrow_left, glyphs.arrow_right
            );
            help_lines.extend(vec![
                Line::from(Span::styled("SEARCH", Style::default().add_modifier(Modifier::BOLD))),
                Line::from(""),
                Line::from("  /pattern       Search forward"),
                Line::from("  ?pattern       Search backward"),
                Line::from("  n              Find next match"),
                Line::from("  N              Find previous match"),
                Line::from(""),
                Line::from("  Searches both sequence names and sequences."),
                Line::from("  Search is case-insensitive."),
                Line::from(""),
                Line::from(Span::styled(tab_hint, Style::default().fg(Color::DarkGray))),
            ]);
        }
        HelpTab::Translation => {
            let ud = format!("{}/{}", glyphs.arrow_up, glyphs.arrow_down);
            let lr = format!("{}/{}", glyphs.arrow_left, glyphs.arrow_right);
            help_lines.extend(vec![
                Line::from(Span::styled("TRANSLATION (NT -> AA)", Style::default().add_modifier(Modifier::BOLD))),
                Line::from(""),
                Line::from("  :asAA          Translate (uses current settings)"),
                Line::from("  :asNT          Switch back to nucleotide view"),
                Line::from("  :setcode       Change genetic code and frame"),
                Line::from(""),
                Line::from("  Default: Standard code (1), frame +1"),
                Line::from("  Settings are remembered between translations."),
                Line::from(""),
                Line::from(Span::styled("SETTINGS DIALOG (:setcode)", Style::default().add_modifier(Modifier::BOLD))),
                Line::from(""),
                Line::from(format!("  {} or j/k     Select genetic code (33 available)", ud)),
                Line::from(format!("  {} or h/l  Select reading frame (+1, +2, +3)", lr)),
                Line::from("  Enter          Confirm       Esc  Cancel"),
            ]);
        }
    }

    let block = Block::default()
        .borders(Borders::ALL)
        .title(" Help ")
        .title_style(Style::default().add_modifier(Modifier::BOLD))
        .style(Style::default().bg(Color::Black));

    let paragraph = Paragraph::new(help_lines).block(block);
    frame.render_widget(paragraph, popup_area);
}

/// Renders the translation settings overlay popup.
fn render_translation_settings_overlay(
    frame: &mut Frame,
    state: &AppState,
    glyphs: &Glyphs,
    area: Rect,
) {
    use crate::genetic_code::GeneticCodes;
    
    // Calculate centered popup dimensions - smaller now since we show only one code
    let popup_width = 56.min(area.width.saturating_sub(4));
    let popup_height = 12.min(area.height.saturating_sub(4));
    
    let popup_x = (area.width.saturating_sub(popup_width)) / 2;
    let popup_y = (area.height.saturating_sub(popup_height)) / 2;
    
    let popup_area = Rect::new(popup_x, popup_y, popup_width, popup_height);

    // Clear the area behind the popup
    frame.render_widget(Clear, popup_area);

    let codes = GeneticCodes::new();
    let all_codes = codes.all();
    let selected_idx = state.translation_settings.selected_code_index;
    let selected_code = all_codes.get(selected_idx);
    
    // Build content
    let mut lines: Vec<Line> = Vec::new();
    
    // Frame selection
    lines.push(Line::from(Span::styled(
        "Reading Frame:",
        Style::default().add_modifier(Modifier::BOLD),
    )));
    
    let frame_spans: Vec<Span> = (0..3)
        .map(|f| {
            let label = format!(" +{} ", f + 1);
            if f == state.translation_settings.selected_frame as usize {
                Span::styled(
                    label,
                    Style::default()
                        .fg(Color::Black)
                        .bg(Color::Cyan)
                        .add_modifier(Modifier::BOLD),
                )
            } else {
                Span::styled(label, Style::default().fg(Color::Gray))
            }
        })
        .collect();
    lines.push(Line::from(vec![
        Span::raw("  "),
        frame_spans[0].clone(),
        Span::raw("  "),
        frame_spans[1].clone(),
        Span::raw("  "),
        frame_spans[2].clone(),
    ]));
    
    lines.push(Line::from(""));
    
    // Genetic code - single line with j/k to change
    lines.push(Line::from(Span::styled(
        "Genetic Code:",
        Style::default().add_modifier(Modifier::BOLD),
    )));
    
    if let Some(code) = selected_code {
        let code_label = format!("  {:2}. {}", code.id, code.name);
        let max_label_len = popup_width as usize;
        let truncated = if code_label.len() > max_label_len.saturating_sub(4) {
            let truncate_len = max_label_len.saturating_sub(7);
            format!("{}...", &code_label[..truncate_len])
        } else {
            code_label
        };
        lines.push(Line::from(Span::styled(
            truncated,
            Style::default()
                .fg(Color::Black)
                .bg(Color::White)
                .add_modifier(Modifier::BOLD),
        )));
    }
    
    // Navigation hint for genetic code
    let nav_hint = format!(
        "  {}/{} to change ({}/{})",
        glyphs.arrow_up,
        glyphs.arrow_down,
        selected_idx + 1,
        all_codes.len()
    );
    lines.push(Line::from(Span::styled(nav_hint, Style::default().fg(Color::DarkGray))));
    
    lines.push(Line::from(""));
    lines.push(Line::from(vec![
        Span::styled(" Enter ", Style::default().fg(Color::Black).bg(Color::DarkGray)),
        Span::styled(" translate  ", Style::default().fg(Color::Gray)),
        Span::styled(" Esc ", Style::default().fg(Color::Black).bg(Color::DarkGray)),
        Span::styled(" cancel ", Style::default().fg(Color::Gray)),
    ]));
    
    lines.push(Line::from(""));
    lines.push(Line::from(Span::styled(
        "  Use :asNT to switch back to nucleotides",
        Style::default().fg(Color::DarkGray),
    )));

    let block = Block::default()
        .borders(Borders::ALL)
        .title(" Translation Settings ")
        .title_style(Style::default().add_modifier(Modifier::BOLD))
        .style(Style::default().bg(Color::Black));

    let paragraph = Paragraph::new(lines).block(block);
    frame.render_widget(paragraph, popup_area);
}

/// Renders a loading overlay with spinner animation.
fn render_loading_overlay(frame: &mut Frame, state: &AppState, area: Rect) {
    use crate::model::LoadingState;
    
    let (message, progress_info) = match &state.loading_state {
        LoadingState::Ready => return, // Nothing to show
        LoadingState::LoadingFile { message, sequences_loaded, .. } => {
            let extra = sequences_loaded
                .map(|n| format!(" ({} sequences)", n))
                .unwrap_or_default();
            (message.clone(), extra)
        }
        LoadingState::Translating { message, sequences_done, total, .. } => {
            let extra = format!(" ({}/{})", sequences_done, total);
            (message.clone(), extra)
        }
    };

    // Calculate centered popup dimensions
    let popup_width = 50.min(area.width.saturating_sub(4));
    let popup_height = 5;
    
    let popup_x = (area.width.saturating_sub(popup_width)) / 2;
    let popup_y = (area.height.saturating_sub(popup_height)) / 2;
    
    let popup_area = Rect::new(popup_x, popup_y, popup_width, popup_height);

    // Clear the area behind the popup
    frame.render_widget(Clear, popup_area);

    // Build content with spinner
    let spinner = state.spinner_char();
    let lines = vec![
        Line::from(""),
        Line::from(vec![
            Span::styled(
                format!("  {} ", spinner),
                Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD),
            ),
            Span::styled(&message, Style::default().fg(Color::White)),
            Span::styled(progress_info, Style::default().fg(Color::DarkGray)),
        ]),
        Line::from(""),
    ];

    let block = Block::default()
        .borders(Borders::ALL)
        .title(" Loading ")
        .title_style(Style::default().add_modifier(Modifier::BOLD))
        .style(Style::default().bg(Color::Black));

    let paragraph = Paragraph::new(lines).block(block);
    frame.render_widget(paragraph, popup_area);
}

/// Renders an error popup overlay.
fn render_error_popup(frame: &mut Frame, error_msg: &str, area: Rect) {
    // Calculate centered popup dimensions
    let popup_width = 60.min(area.width.saturating_sub(4));
    let popup_height = 7;
    
    let popup_x = (area.width.saturating_sub(popup_width)) / 2;
    let popup_y = (area.height.saturating_sub(popup_height)) / 2;
    
    let popup_area = Rect::new(popup_x, popup_y, popup_width, popup_height);

    // Clear the area behind the popup
    frame.render_widget(Clear, popup_area);

    // Wrap error message if needed
    let max_line_len = (popup_width - 4) as usize;
    let wrapped_msg = textwrap::wrap(error_msg, max_line_len);
    
    let mut lines = vec![Line::from("")];
    for line in wrapped_msg.iter().take(3) {
        lines.push(Line::from(Span::styled(
            format!("  {}", line),
            Style::default().fg(Color::Red),
        )));
    }
    lines.push(Line::from(""));
    lines.push(Line::from(Span::styled(
        "  Press any key to continue...",
        Style::default().fg(Color::DarkGray),
    )));

    let block = Block::default()
        .borders(Borders::ALL)
        .border_style(Style::default().fg(Color::Red))
        .title(" Error ")
        .title_style(Style::default().fg(Color::Red).add_modifier(Modifier::BOLD))
        .style(Style::default().bg(Color::Black));

    let paragraph = Paragraph::new(lines).block(block);
    frame.render_widget(paragraph, popup_area);
}

/// Renders the file browser overlay.
fn render_file_browser(
    frame: &mut Frame,
    browser: &crate::model::FileBrowserState,
    glyphs: &Glyphs,
    area: Rect,
) {
    // Calculate centered popup dimensions
    let popup_width = 70.min(area.width.saturating_sub(4));
    let popup_height = 20.min(area.height.saturating_sub(4));
    
    let popup_x = (area.width.saturating_sub(popup_width)) / 2;
    let popup_y = (area.height.saturating_sub(popup_height)) / 2;
    
    let popup_area = Rect::new(popup_x, popup_y, popup_width, popup_height);

    // Clear the area behind the popup
    frame.render_widget(Clear, popup_area);

    // Calculate visible height for entries (excluding borders and title)
    let visible_height = (popup_height.saturating_sub(4)) as usize;

    // Build entry lines
    let mut lines: Vec<Line> = Vec::new();
    
    // Show current directory path as first line
    let dir_str = browser.current_dir.display().to_string();
    let max_path_len = (popup_width - 4) as usize;
    let display_path = if dir_str.len() > max_path_len {
        format!("...{}", &dir_str[dir_str.len() - max_path_len + 3..])
    } else {
        dir_str
    };
    lines.push(Line::from(Span::styled(
        display_path,
        Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD),
    )));
    lines.push(Line::from(glyphs.h_separator.repeat((popup_width - 2) as usize)));

    // Determine scroll window
    let start_idx = browser.scroll_offset;
    let end_idx = (start_idx + visible_height.saturating_sub(2)).min(browser.entries.len());

    for idx in start_idx..end_idx {
        let entry = &browser.entries[idx];
        let is_selected = idx == browser.selected;

        let (prefix, style) = if entry.is_dir {
            (glyphs.dir_prefix, Style::default().fg(Color::Yellow))
        } else {
            (glyphs.file_prefix, Style::default().fg(Color::White))
        };

        let name_style = if is_selected {
            style.bg(Color::DarkGray).add_modifier(Modifier::BOLD)
        } else {
            style
        };

        // Truncate name if too long
        let max_name_len = (popup_width - 6) as usize;
        let display_name = if entry.name.len() > max_name_len {
            let truncate_len = max_name_len.saturating_sub(3);
            format!("{}...", &entry.name[..truncate_len])
        } else {
            entry.name.clone()
        };

        let line_content = format!("{}{}", prefix, display_name);
        
        // Pad to full width if selected (for highlight)
        let padded = if is_selected {
            format!("{:<width$}", line_content, width = (popup_width - 2) as usize)
        } else {
            line_content
        };

        lines.push(Line::from(Span::styled(padded, name_style)));
    }

    // Show hint at bottom
    lines.push(Line::from(""));
    let nav_hint = format!(
        " {}/{}:Navigate  Enter:Select  Backspace:Parent  Esc:Quit",
        glyphs.arrow_up, glyphs.arrow_down
    );
    lines.push(Line::from(Span::styled(
        nav_hint,
        Style::default().fg(Color::DarkGray),
    )));

    let title = format!(" Select File - {} ", browser.error_message);
    let block = Block::default()
        .borders(Borders::ALL)
        .title(title)
        .title_style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD))
        .style(Style::default().bg(Color::Black));

    let paragraph = Paragraph::new(lines).block(block);
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
