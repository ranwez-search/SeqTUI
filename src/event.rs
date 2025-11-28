//! Keyboard event handling.
//!
//! This module manages keyboard input with Vim-style navigation:
//! - `h`: move left
//! - `j`: move down
//! - `k`: move up
//! - `l`: move right
//! - `0` or `Home`: go to first column
//! - `$` or `End`: go to last column
//! - `g0`: go to first visible column
//! - `gm`: go to middle of visible area
//! - `g$`: go to last visible column
//! - `<number>|`: go to column (e.g., `50|` goes to column 50)
//! - `:`: enter command mode
//!   - `:q` or `:quit`: quit the application
//!   - `:h` or `:help`: show help
//!   - `:<number>`: go to column
//! - `/`: search forward
//! - `?`: search backward
//! - `n`: find next match
//! - `N`: find previous match

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
    /// Enter search mode (forward with /)
    EnterSearchMode,
    /// Enter search mode (backward with ?)
    EnterSearchBackward,
    /// Add character to search buffer
    SearchChar(char),
    /// Execute search
    ExecuteSearch,
    /// Cancel search mode
    CancelSearch,
    /// Backspace in search mode
    SearchBackspace,
    /// Find next match (n)
    FindNext,
    /// Find previous match (N)
    FindPrevious,
    /// Dismiss the help overlay
    DismissHelp,
    /// Go to first column (0 or Home)
    GotoFirstColumn,
    /// Go to last column ($ or End)
    GotoLastColumn,
    /// Go to first visible column (g0)
    GotoFirstVisibleColumn,
    /// Go to middle of visible area (gm)
    GotoMiddleVisibleColumn,
    /// Go to last visible column (g$)
    GotoLastVisibleColumn,
    /// Go to specific column (number|)
    GotoColumn(usize),
    /// Pending 'g' key for g-commands
    PendingG,
    /// Accumulate a digit for number prefix
    AccumulateDigit(char),
    /// Execute pending number with | (go to column)
    ExecuteGotoColumn,
    /// Translation settings: move selection up
    TranslationUp,
    /// Translation settings: move selection down
    TranslationDown,
    /// Translation settings: select previous frame
    TranslationFrameLeft,
    /// Translation settings: select next frame
    TranslationFrameRight,
    /// Translation settings: confirm and translate
    TranslationConfirm,
    /// Translation settings: cancel
    TranslationCancel,
    /// Move half page up (Ctrl+U)
    HalfPageUp,
    /// Move half page down (Ctrl+D)
    HalfPageDown,
    /// Move full page up (PageUp)
    PageUp,
    /// Move full page down (PageDown)
    PageDown,
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
pub fn handle_event(event: Event, mode: &AppMode, show_help: bool, pending_g: bool, has_number_prefix: bool) -> Action {
    match event {
        Event::Key(key_event) => handle_key_event(key_event, mode, show_help, pending_g, has_number_prefix),
        Event::Resize(width, height) => Action::Resize(width, height),
        _ => Action::None,
    }
}

/// Handles a key event based on the current application mode.
fn handle_key_event(key: KeyEvent, mode: &AppMode, show_help: bool, pending_g: bool, has_number_prefix: bool) -> Action {
    // If help is shown, any key dismisses it
    if show_help {
        return Action::DismissHelp;
    }

    // Handle pending g-command
    if pending_g {
        return handle_g_command(key);
    }

    match mode {
        AppMode::Normal => handle_normal_mode(key, has_number_prefix),
        AppMode::Command(_) => handle_command_mode(key),
        AppMode::Search(_) | AppMode::SearchBackward(_) => handle_search_mode(key),
        AppMode::TranslationSettings => handle_translation_settings_mode(key),
    }
}

/// Handles g-prefix commands (g0, gm, g$).
fn handle_g_command(key: KeyEvent) -> Action {
    match key.code {
        KeyCode::Char('0') => Action::GotoFirstVisibleColumn,
        KeyCode::Char('m') => Action::GotoMiddleVisibleColumn,
        KeyCode::Char('$') => Action::GotoLastVisibleColumn,
        _ => Action::None, // Unknown g-command, cancel
    }
}

/// Handles key events in normal mode (Vim-style navigation).
fn handle_normal_mode(key: KeyEvent, has_number_prefix: bool) -> Action {
    // Handle Ctrl+C for emergency quit
    if key.modifiers.contains(KeyModifiers::CONTROL) && key.code == KeyCode::Char('c') {
        return Action::Quit;
    }

    // Handle Ctrl+U for half page up
    if key.modifiers.contains(KeyModifiers::CONTROL) && key.code == KeyCode::Char('u') {
        return Action::HalfPageUp;
    }

    // Handle Ctrl+D for half page down
    if key.modifiers.contains(KeyModifiers::CONTROL) && key.code == KeyCode::Char('d') {
        return Action::HalfPageDown;
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

        // Jump to start/end of line (Vim-style)
        // '0' goes to first column only if not accumulating a number (e.g., 500|)
        KeyCode::Char('0') if !has_number_prefix => Action::GotoFirstColumn,
        KeyCode::Char('0') => Action::AccumulateDigit('0'),
        KeyCode::Char('$') => Action::GotoLastColumn,
        KeyCode::Home => Action::GotoFirstColumn,
        KeyCode::End => Action::GotoLastColumn,

        // Page navigation
        KeyCode::PageUp => Action::PageUp,
        KeyCode::PageDown => Action::PageDown,

        // g-prefix commands (handled via pending state)
        KeyCode::Char('g') => Action::PendingG,

        // Number prefix for <number>| command
        KeyCode::Char(c @ '1'..='9') => Action::AccumulateDigit(c),
        KeyCode::Char('|') => Action::ExecuteGotoColumn,

        // Command mode
        KeyCode::Char(':') => Action::EnterCommandMode,

        // Search mode (Vim-style)
        KeyCode::Char('/') => Action::EnterSearchMode,
        KeyCode::Char('?') => Action::EnterSearchBackward,
        KeyCode::Char('n') => Action::FindNext,
        KeyCode::Char('N') => Action::FindPrevious,

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

/// Handles key events in search mode.
fn handle_search_mode(key: KeyEvent) -> Action {
    match key.code {
        KeyCode::Enter => Action::ExecuteSearch,
        KeyCode::Esc => Action::CancelSearch,
        KeyCode::Backspace => Action::SearchBackspace,
        KeyCode::Char(c) => Action::SearchChar(c),
        _ => Action::None,
    }
}

/// Handles key events in translation settings mode.
fn handle_translation_settings_mode(key: KeyEvent) -> Action {
    match key.code {
        // Navigation
        KeyCode::Char('j') | KeyCode::Down => Action::TranslationDown,
        KeyCode::Char('k') | KeyCode::Up => Action::TranslationUp,
        KeyCode::Char('h') | KeyCode::Left => Action::TranslationFrameLeft,
        KeyCode::Char('l') | KeyCode::Right => Action::TranslationFrameRight,
        // Frame selection with number keys
        KeyCode::Char('1') => Action::TranslationFrameLeft, // Will cycle, or we can set directly
        KeyCode::Char('2') => Action::TranslationFrameRight,
        KeyCode::Char('3') => Action::TranslationFrameRight,
        // Confirm / Cancel
        KeyCode::Enter => Action::TranslationConfirm,
        KeyCode::Esc | KeyCode::Char('q') => Action::TranslationCancel,
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
        Action::EnterSearchMode => {
            state.enter_search_mode(false);
        }
        Action::EnterSearchBackward => {
            state.enter_search_mode(true);
        }
        Action::SearchChar(c) => {
            state.search_input(c);
        }
        Action::ExecuteSearch => {
            state.execute_search();
        }
        Action::CancelSearch => {
            state.cancel_search();
        }
        Action::SearchBackspace => {
            state.search_backspace();
        }
        Action::FindNext => {
            state.find_next();
        }
        Action::FindPrevious => {
            state.find_previous();
        }
        Action::DismissHelp => {
            state.dismiss_help();
        }
        Action::GotoFirstColumn => {
            state.goto_first_column();
        }
        Action::GotoLastColumn => {
            state.goto_last_column();
        }
        Action::GotoFirstVisibleColumn => {
            state.goto_first_visible_column();
        }
        Action::GotoMiddleVisibleColumn => {
            state.goto_middle_visible_column();
        }
        Action::GotoLastVisibleColumn => {
            state.goto_last_visible_column();
        }
        Action::GotoColumn(col) => {
            state.goto_column(col);
        }
        Action::PendingG => {
            state.set_pending_g();
        }
        Action::AccumulateDigit(c) => {
            state.accumulate_digit(c);
        }
        Action::ExecuteGotoColumn => {
            state.execute_goto_column();
        }
        Action::TranslationUp => {
            state.translation_settings_up();
        }
        Action::TranslationDown => {
            state.translation_settings_down();
        }
        Action::TranslationFrameLeft => {
            state.translation_settings_frame_left();
        }
        Action::TranslationFrameRight => {
            state.translation_settings_frame_right();
        }
        Action::TranslationConfirm => {
            state.confirm_translation_settings();
        }
        Action::TranslationCancel => {
            state.cancel_translation_settings();
        }
        Action::HalfPageUp => {
            state.half_page_up();
        }
        Action::HalfPageDown => {
            state.half_page_down();
        }
        Action::PageUp => {
            state.page_up();
        }
        Action::PageDown => {
            state.page_down();
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
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::MoveLeft);

        let key = KeyEvent::new(KeyCode::Char('j'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::MoveDown);

        let key = KeyEvent::new(KeyCode::Char('k'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::MoveUp);

        let key = KeyEvent::new(KeyCode::Char('l'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::MoveRight);
    }

    #[test]
    fn test_enter_command_mode() {
        let mode = AppMode::Normal;
        let key = KeyEvent::new(KeyCode::Char(':'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::EnterCommandMode);
    }

    #[test]
    fn test_command_mode_input() {
        let mode = AppMode::Command(String::new());

        let key = KeyEvent::new(KeyCode::Char('q'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::CommandChar('q'));

        let key = KeyEvent::new(KeyCode::Enter, KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::ExecuteCommand);

        let key = KeyEvent::new(KeyCode::Esc, KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::CancelCommand);
    }

    #[test]
    fn test_ctrl_c_quit() {
        let mode = AppMode::Normal;
        let key = KeyEvent::new(KeyCode::Char('c'), KeyModifiers::CONTROL);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::Quit);
    }

    #[test]
    fn test_search_mode_keys() {
        let mode = AppMode::Normal;

        // Test entering search modes
        let key = KeyEvent::new(KeyCode::Char('/'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::EnterSearchMode);

        let key = KeyEvent::new(KeyCode::Char('?'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::EnterSearchBackward);

        // Test find next/previous
        let key = KeyEvent::new(KeyCode::Char('n'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::FindNext);

        let key = KeyEvent::new(KeyCode::Char('N'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::FindPrevious);
    }

    #[test]
    fn test_search_mode_input() {
        let mode = AppMode::Search(String::new());

        let key = KeyEvent::new(KeyCode::Char('A'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::SearchChar('A'));

        let key = KeyEvent::new(KeyCode::Enter, KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::ExecuteSearch);

        let key = KeyEvent::new(KeyCode::Esc, KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::CancelSearch);

        let key = KeyEvent::new(KeyCode::Backspace, KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::SearchBackspace);
    }

    #[test]
    fn test_dismiss_help() {
        let mode = AppMode::Normal;
        // Any key when help is shown should dismiss help
        let key = KeyEvent::new(KeyCode::Char('x'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, true, false, false), Action::DismissHelp);

        let key = KeyEvent::new(KeyCode::Esc, KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, true, false, false), Action::DismissHelp);
    }

    #[test]
    fn test_jump_navigation() {
        let mode = AppMode::Normal;

        // Test 0 key (without number prefix -> go to first column)
        let key = KeyEvent::new(KeyCode::Char('0'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::GotoFirstColumn);

        // Test 0 key WITH number prefix (e.g., typing 500|) -> accumulate digit
        let key = KeyEvent::new(KeyCode::Char('0'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, true), Action::AccumulateDigit('0'));

        let key = KeyEvent::new(KeyCode::Char('$'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::GotoLastColumn);

        // Test Home and End keys
        let key = KeyEvent::new(KeyCode::Home, KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::GotoFirstColumn);

        let key = KeyEvent::new(KeyCode::End, KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::GotoLastColumn);
    }

    #[test]
    fn test_g_commands() {
        let mode = AppMode::Normal;

        // Test g prefix
        let key = KeyEvent::new(KeyCode::Char('g'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::PendingG);

        // Test g0, gm, g$ (with pending_g = true)
        let key = KeyEvent::new(KeyCode::Char('0'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, true, false), Action::GotoFirstVisibleColumn);

        let key = KeyEvent::new(KeyCode::Char('m'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, true, false), Action::GotoMiddleVisibleColumn);

        let key = KeyEvent::new(KeyCode::Char('$'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, true, false), Action::GotoLastVisibleColumn);

        // Unknown g-command returns None
        let key = KeyEvent::new(KeyCode::Char('x'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, true, false), Action::None);
    }

    #[test]
    fn test_number_prefix() {
        let mode = AppMode::Normal;

        // Test digit accumulation
        let key = KeyEvent::new(KeyCode::Char('5'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::AccumulateDigit('5'));

        // Test | to execute goto column
        let key = KeyEvent::new(KeyCode::Char('|'), KeyModifiers::NONE);
        assert_eq!(handle_key_event(key, &mode, false, false, false), Action::ExecuteGotoColumn);
    }
}
