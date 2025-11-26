//! Keyboard event handling.
//!
//! This module manages keyboard input with Vim-style navigation:
//! - `h`: move left
//! - `j`: move down
//! - `k`: move up
//! - `l`: move right
//! - `:`: enter command mode
//! - `q` in command mode: quit
//!
//! The design supports future extensions like:
//! - `/` for pattern search
//! - `g` for go-to commands
//! - Number prefixes for repeated movements

use crossterm::event::{self, Event, KeyCode, KeyEvent, KeyModifiers};
use std::time::Duration;

use crate::model::{AppMode, AppState};

/// Actions that can be triggered by keyboard input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Action {
    /// No action (key not recognized)
    None,
    /// Quit the application
    Quit,
    /// Move cursor up
    MoveUp,
    /// Move cursor down
    MoveDown,
    /// Move cursor left
    MoveLeft,
    /// Move cursor right
    MoveRight,
    /// Enter command mode
    EnterCommandMode,
    /// Add character to command buffer
    CommandChar(char),
    /// Execute current command
    ExecuteCommand,
    /// Cancel command mode
    CancelCommand,
    /// Backspace in command mode
    CommandBackspace,
    /// Resize event (terminal resized)
    Resize(u16, u16),
}

/// Polls for keyboard events with a timeout.
///
/// Returns `None` if no event occurred within the timeout.
pub fn poll_event(timeout: Duration) -> Option<Event> {
    if event::poll(timeout).ok()? {
        event::read().ok()
    } else {
        None
    }
}

/// Converts a crossterm event to an Action based on current app mode.
pub fn handle_event(event: Event, mode: &AppMode) -> Action {
    match event {
        Event::Key(key_event) => handle_key_event(key_event, mode),
        Event::Resize(width, height) => Action::Resize(width, height),
        _ => Action::None,
    }
}

/// Handles a key event based on the current application mode.
fn handle_key_event(key: KeyEvent, mode: &AppMode) -> Action {
    match mode {
        AppMode::Normal => handle_normal_mode(key),
        AppMode::Command(_) => handle_command_mode(key),
    }
}

/// Handles key events in normal mode (Vim-style navigation).
fn handle_normal_mode(key: KeyEvent) -> Action {
    // Handle Ctrl+C for emergency quit
    if key.modifiers.contains(KeyModifiers::CONTROL) && key.code == KeyCode::Char('c') {
        return Action::Quit;
    }

    match key.code {
        // Vim-style navigation (as per spec: j=down, k=up, l=right, h=left)
        KeyCode::Char('j') => Action::MoveDown,
        KeyCode::Char('k') => Action::MoveUp,
        KeyCode::Char('l') => Action::MoveRight,
        KeyCode::Char('h') => Action::MoveLeft,

        // Alternative arrow keys for convenience
        KeyCode::Up => Action::MoveUp,
        KeyCode::Down => Action::MoveDown,
        KeyCode::Right => Action::MoveRight,
        KeyCode::Left => Action::MoveLeft,

        // Command mode
        KeyCode::Char(':') => Action::EnterCommandMode,

        // Quick quit with 'q' in normal mode (optional convenience)
        // Commented out to strictly follow spec (use :q)
        // KeyCode::Char('q') => Action::Quit,

        _ => Action::None,
    }
}

/// Handles key events in command mode.
fn handle_command_mode(key: KeyEvent) -> Action {
    match key.code {
        KeyCode::Enter => Action::ExecuteCommand,
        KeyCode::Esc => Action::CancelCommand,
        KeyCode::Backspace => Action::CommandBackspace,
        KeyCode::Char(c) => Action::CommandChar(c),
        _ => Action::None,
    }
}

/// Applies an action to the application state.
///
/// Returns `true` if the application should continue, `false` if it should quit.
pub fn apply_action(state: &mut AppState, action: Action) -> bool {
    match action {
        Action::None => {}
        Action::Quit => {
            state.should_quit = true;
        }
        Action::MoveUp => {
            state.move_up();
        }
        Action::MoveDown => {
            state.move_down();
        }
        Action::MoveLeft => {
            state.move_left();
        }
        Action::MoveRight => {
            state.move_right();
        }
        Action::EnterCommandMode => {
            state.enter_command_mode();
        }
        Action::CommandChar(c) => {
            state.command_input(c);
        }
        Action::ExecuteCommand => {
            state.execute_command();
        }
        Action::CancelCommand => {
            state.cancel_command();
        }
        Action::CommandBackspace => {
            state.command_backspace();
        }
        Action::Resize(_, _) => {
            // Resize is handled in the main loop with actual terminal dimensions
        }
    }

    !state.should_quit
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normal_mode_navigation() {
        let mode = AppMode::Normal;

        // Test movement keys (Vim-style: h=left, j=down, k=up, l=right)
        let key = KeyEvent::new(KeyCode::Char('h'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode), Action::MoveLeft);

        let key = KeyEvent::new(KeyCode::Char('j'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode), Action::MoveDown);

        let key = KeyEvent::new(KeyCode::Char('k'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode), Action::MoveUp);

        let key = KeyEvent::new(KeyCode::Char('l'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode), Action::MoveRight);
    }

    #[test]
    fn test_enter_command_mode() {
        let mode = AppMode::Normal;
        let key = KeyEvent::new(KeyCode::Char(':'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode), Action::EnterCommandMode);
    }

    #[test]
    fn test_command_mode_input() {
        let mode = AppMode::Command(String::new());

        let key = KeyEvent::new(KeyCode::Char('q'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode), Action::CommandChar('q'));

        let key = KeyEvent::new(KeyCode::Enter, KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode), Action::ExecuteCommand);

        let key = KeyEvent::new(KeyCode::Esc, KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode), Action::CancelCommand);
    }

    #[test]
    fn test_ctrl_c_quit() {
        let mode = AppMode::Normal;
        let key = KeyEvent::new(KeyCode::Char('c'), KeyModifiers::CONTROL);
        assert_eq!(handle_key_event(key, &mode), Action::Quit);
    }
}
