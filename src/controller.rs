//! Application controller.
//!
//! This module orchestrates the main application loop:
//! - Terminal initialization and cleanup
//! - Event polling and handling
//! - State updates and rendering
//!
//! The design supports future extensions like:
//! - Multiple views (file browser, help screen)
//! - Async operations (background file loading)
//! - Plugin system

use std::io::{self, Stdout};
use std::time::Duration;

use anyhow::Result;
use crossterm::{
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{backend::CrosstermBackend, Terminal};

use crate::event::{apply_action, handle_event, poll_event};
use crate::model::AppState;
use crate::ui::{calculate_visible_dimensions, render};

/// The main application controller.
pub struct App {
    /// Terminal backend
    terminal: Terminal<CrosstermBackend<Stdout>>,
    /// Application state
    state: AppState,
    /// Event poll timeout
    tick_rate: Duration,
}

impl App {
    /// Creates a new application with the given state.
    pub fn new(state: AppState) -> Result<Self> {
        // Setup terminal
        enable_raw_mode()?;
        let mut stdout = io::stdout();
        execute!(stdout, EnterAlternateScreen)?;
        let backend = CrosstermBackend::new(stdout);
        let terminal = Terminal::new(backend)?;

        Ok(Self {
            terminal,
            state,
            tick_rate: Duration::from_millis(50),
        })
    }

    /// Runs the main application loop.
    pub fn run(&mut self) -> Result<()> {
        // Initial viewport setup
        self.update_viewport_size()?;

        loop {
            // Render
            self.terminal.draw(|frame| {
                render(frame, &self.state);
            })?;

            // Handle events
            if let Some(event) = poll_event(self.tick_rate) {
                let has_number_prefix = !self.state.number_buffer.is_empty();
                let action = handle_event(event, &self.state.mode, self.state.show_help, self.state.pending_g, self.state.pending_z, has_number_prefix);

                // Handle resize specially to update viewport
                if let crate::event::Action::Resize(_, _) = action {
                    self.update_viewport_size()?;
                }

                apply_action(&mut self.state, action);

                if self.state.should_quit {
                    break;
                }
            }
        }

        Ok(())
    }

    /// Updates the viewport size based on terminal dimensions.
    fn update_viewport_size(&mut self) -> Result<()> {
        let size = self.terminal.size()?;
        let (visible_rows, visible_cols) = calculate_visible_dimensions(size.width, size.height);
        self.state.update_viewport_size(visible_rows, visible_cols);
        Ok(())
    }
}

impl Drop for App {
    fn drop(&mut self) {
        // Restore terminal
        let _ = disable_raw_mode();
        let _ = execute!(self.terminal.backend_mut(), LeaveAlternateScreen);
        let _ = self.terminal.show_cursor();
    }
}

/// Convenience function to run the application with an alignment file.
pub fn run_app(state: AppState) -> Result<()> {
    let mut app = App::new(state)?;
    app.run()
}

#[cfg(test)]
mod tests {
    use crate::model::{Alignment, Sequence};

    use super::*;

    #[test]
    fn test_app_state_creation() {
        let seqs = vec![
            Sequence::new("seq1", "ACGT"),
            Sequence::new("seq2", "TGCA"),
        ];
        let alignment = Alignment::new(seqs);
        let state = AppState::new(alignment, "test".to_string());

        assert_eq!(state.alignment.sequence_count(), 2);
        assert!(!state.should_quit);
    }
}
