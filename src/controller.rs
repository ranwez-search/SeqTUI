//! Application controller.
//!
//! This module orchestrates the main application loop:
//! - Terminal initialization and cleanup
//! - Event polling and handling
//! - State updates and rendering
//! - Background file loading with progress updates
//! - Background translation with progress updates
//!
//! The design supports:
//! - Async file loading (TUI shows spinner while parsing)
//! - Async translation (TUI shows progress while translating)
//! - Multiple views (help screen, translation settings)

use std::io::{self, Stdout};
use std::path::PathBuf;
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread;
use std::time::Duration;

use anyhow::Result;
use crossterm::{
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{backend::CrosstermBackend, Terminal};

use crate::event::{apply_action, handle_event, poll_event, ActionResult};
use crate::formats::{parse_file_with_options, FileFormat};
use crate::genetic_code::GeneticCodes;
use crate::model::{Alignment, AppState, LoadingState, Sequence, SequenceType};
use crate::ui::{calculate_visible_dimensions, render};

/// Messages sent from the background loading thread.
pub enum LoadMessage {
    /// Loading completed successfully
    Complete(Alignment),
    /// Loading failed with an error (includes file path for context)
    Error { message: String, path: PathBuf },
    /// Progress update (for future use with streaming parsers)
    Progress { sequences_loaded: usize },
}

/// Messages sent from the background translation thread.
pub enum TranslateMessage {
    /// Translation progress update
    Progress { sequences_done: usize, total: usize },
    /// Translation completed successfully
    Complete(Alignment),
}

/// The main application controller.
pub struct App {
    /// Terminal backend
    terminal: Terminal<CrosstermBackend<Stdout>>,
    /// Application state
    state: AppState,
    /// Event poll timeout
    tick_rate: Duration,
    /// Receiver for background loading messages
    load_receiver: Option<Receiver<LoadMessage>>,
    /// Receiver for background translation messages
    translate_receiver: Option<Receiver<TranslateMessage>>,
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
            load_receiver: None,
            translate_receiver: None,
        })
    }

    /// Creates a new application and starts loading a file in the background.
    /// Optional `preset_translation` is (genetic_code_id, reading_frame) to preset translation settings.
    pub fn new_with_background_load(
        file_path: PathBuf,
        forced_format: Option<FileFormat>,
        preset_translation: Option<(u8, u8)>,
        fancy_ui: bool,
    ) -> Result<Self> {
        // Extract file name for display
        let file_name = file_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("alignment")
            .to_string();

        // Create initial state in loading mode
        let mut state = AppState::new_loading(file_name, file_path.clone());
        state.fancy_ui = fancy_ui;
        
        // Apply preset translation settings if provided
        if let Some((genetic_code, reading_frame)) = preset_translation {
            // Find the index of the genetic code in the list
            let codes = GeneticCodes::new();
            let code_index = codes.all().iter().position(|c| c.id == genetic_code).unwrap_or(0);
            
            state.translation_settings.genetic_code_id = genetic_code;
            state.translation_settings.selected_code_index = code_index;
            // Convert from 1-based to 0-based frame
            let frame = (reading_frame.saturating_sub(1) as usize).min(2);
            state.translation_settings.frame = frame;
            state.translation_settings.selected_frame = frame;
            state.translation_settings.has_translated = true; // Skip the settings dialog
        }

        // Setup terminal
        enable_raw_mode()?;
        let mut stdout = io::stdout();
        execute!(stdout, EnterAlternateScreen)?;
        let backend = CrosstermBackend::new(stdout);
        let terminal = Terminal::new(backend)?;

        // Create channel for background loading
        let (tx, rx): (Sender<LoadMessage>, Receiver<LoadMessage>) = mpsc::channel();

        // Spawn background thread for loading
        let path_for_error = file_path.clone();
        thread::spawn(move || {
            match parse_file_with_options(&file_path, forced_format) {
                Ok(alignment) => {
                    let _ = tx.send(LoadMessage::Complete(alignment));
                }
                Err(e) => {
                    let _ = tx.send(LoadMessage::Error { 
                        message: e.to_string(), 
                        path: path_for_error,
                    });
                }
            }
        });

        Ok(Self {
            terminal,
            state,
            tick_rate: Duration::from_millis(50),
            load_receiver: Some(rx),
            translate_receiver: None,
        })
    }

    /// Starts background translation of the current alignment.
    fn start_background_translation(&mut self) {
        // Get translation parameters
        let genetic_code_id = self.state.translation_settings.genetic_code_id;
        let frame = self.state.translation_settings.frame;
        
        // Get the genetic code before spawning (avoids recreating all codes in thread)
        let codes = GeneticCodes::new();
        let code = codes.get(genetic_code_id)
            .unwrap_or_else(|| codes.default_code())
            .clone(); // Clone just the one code we need (64 bytes + name)
        
        // Clone the sequences for the background thread
        // Unfortunately we need to copy the data since we can't share references across threads
        let sequences: Vec<(String, Vec<u8>)> = self.state.alignment.sequences
            .iter()
            .map(|seq| (seq.id.clone(), seq.clone_bytes()))
            .collect();
        
        let total = sequences.len();
        
        // Set loading state
        self.state.loading_state = LoadingState::Translating {
            message: "Translating...".to_string(),
            sequences_done: 0,
            total,
        };
        
        // Create channel for translation messages
        let (tx, rx): (Sender<TranslateMessage>, Receiver<TranslateMessage>) = mpsc::channel();
        self.translate_receiver = Some(rx);
        
        // Spawn background thread for translation
        thread::spawn(move || {
            let mut translated_seqs = Vec::with_capacity(total);
            let progress_interval = (total / 20).max(1); // Update every 5%
            
            for (i, (id, data)) in sequences.into_iter().enumerate() {
                let aa_data = code.translate_sequence(&data, frame);
                translated_seqs.push(Sequence::from_bytes(id, aa_data));
                
                // Send progress update periodically
                if (i + 1) % progress_interval == 0 || i + 1 == total {
                    let _ = tx.send(TranslateMessage::Progress {
                        sequences_done: i + 1,
                        total,
                    });
                }
            }
            
            translated_seqs.shrink_to_fit();
            let mut alignment = Alignment::new(translated_seqs);
            alignment.sequence_type = SequenceType::AMINO_ACID;
            let _ = tx.send(TranslateMessage::Complete(alignment));
        });
    }

    /// Runs the main application loop.
    pub fn run(&mut self) -> Result<()> {
        // Initial viewport setup
        self.update_viewport_size()?;

        loop {
            // Check for background loading messages (non-blocking)
            if let Some(ref rx) = self.load_receiver {
                match rx.try_recv() {
                    Ok(LoadMessage::Complete(alignment)) => {
                        self.state.set_alignment(alignment);
                        self.load_receiver = None; // Done loading
                    }
                    Ok(LoadMessage::Error { message, path }) => {
                        self.state.set_loading_error(message, Some(path));
                        self.load_receiver = None;
                    }
                    Ok(LoadMessage::Progress { sequences_loaded }) => {
                        // Update progress in loading state
                        if let LoadingState::LoadingFile { sequences_loaded: ref mut sl, .. } = self.state.loading_state {
                            *sl = Some(sequences_loaded);
                        }
                    }
                    Err(mpsc::TryRecvError::Empty) => {
                        // No message yet, continue
                    }
                    Err(mpsc::TryRecvError::Disconnected) => {
                        // Sender dropped without sending - should not happen
                        self.state.set_loading_error("Loading thread terminated unexpectedly".to_string(), None);
                        self.load_receiver = None;
                    }
                }
            }

            // Check for background translation messages (non-blocking)
            if let Some(ref rx) = self.translate_receiver {
                match rx.try_recv() {
                    Ok(TranslateMessage::Progress { sequences_done, total }) => {
                        if let LoadingState::Translating { sequences_done: ref mut sd, total: ref mut t, .. } = self.state.loading_state {
                            *sd = sequences_done;
                            *t = total;
                        }
                    }
                    Ok(TranslateMessage::Complete(alignment)) => {
                        self.state.set_translated_alignment(alignment);
                        self.translate_receiver = None;
                    }
                    Err(mpsc::TryRecvError::Empty) => {
                        // No message yet, continue
                    }
                    Err(mpsc::TryRecvError::Disconnected) => {
                        // Sender dropped without sending - should not happen
                        self.state.loading_state = LoadingState::Ready;
                        self.state.status_message = Some("Translation failed unexpectedly".to_string());
                        self.translate_receiver = None;
                    }
                }
            }

            // Tick spinner animation if loading
            if self.state.loading_state.is_loading() {
                self.state.tick_spinner();
            }

            // Render
            self.terminal.draw(|frame| {
                render(frame, &self.state);
            })?;

            // Handle events
            if let Some(event) = poll_event(self.tick_rate) {
                let has_number_prefix = !self.state.number_buffer.is_empty();
                let has_error_popup = self.state.error_popup.is_some();
                let has_file_browser = self.state.file_browser.is_some();
                let action = handle_event(
                    event, 
                    &self.state.mode, 
                    self.state.show_help, 
                    self.state.pending_g, 
                    self.state.pending_z, 
                    has_number_prefix,
                    has_error_popup,
                    has_file_browser,
                );

                // Handle resize specially to update viewport
                if let crate::event::Action::Resize(_, _) = action {
                    self.update_viewport_size()?;
                }

                let result = apply_action(&mut self.state, action);

                // Check if translation should be started
                match result {
                    ActionResult::StartTranslation => {
                        self.start_background_translation();
                    }
                    ActionResult::LoadFile(path) => {
                        self.start_background_load(path);
                    }
                    ActionResult::Continue => {}
                }

                if self.state.should_quit {
                    break;
                }
            }
        }

        Ok(())
    }

    /// Starts a new background load for the given file path.
    fn start_background_load(&mut self, file_path: PathBuf) {
        // Close file browser
        self.state.close_file_browser();
        
        // Extract file name for display
        let file_name = file_path.file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("Unknown")
            .to_string();
        
        // Reset state for loading
        self.state.loading_state = LoadingState::LoadingFile {
            path: file_path.clone(),
            message: format!("Loading {}...", file_name),
            sequences_loaded: None,
        };
        self.state.status_message = None;
        
        // Create channel for background loading
        let (tx, rx): (Sender<LoadMessage>, Receiver<LoadMessage>) = mpsc::channel();
        
        // Store the path for error handling
        let path_for_error = file_path.clone();
        
        // Spawn background thread for loading
        thread::spawn(move || {
            match parse_file_with_options(&file_path, None) {
                Ok(alignment) => {
                    let _ = tx.send(LoadMessage::Complete(alignment));
                }
                Err(e) => {
                    let _ = tx.send(LoadMessage::Error { 
                        message: e.to_string(), 
                        path: path_for_error,
                    });
                }
            }
        });
        
        self.load_receiver = Some(rx);
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

/// Convenience function to run the application with the file browser open.
/// Used when no file is provided on command line.
pub fn run_app_with_file_browser(fancy_ui: bool) -> Result<()> {
    use crate::model::{Alignment, FileBrowserState};
    
    let start_dir = std::env::current_dir().unwrap_or_else(|_| std::path::PathBuf::from("."));
    let mut state = AppState::new(Alignment::new(vec![]), "No file".to_string());
    state.fancy_ui = fancy_ui;
    state.file_browser = Some(FileBrowserState::new(start_dir, "Select a sequence file".to_string()));
    
    let mut app = App::new(state)?;
    app.run()
}

/// Convenience function to run the application with the file browser open in a given directory.
pub fn run_app_with_file_browser_at(start_dir: PathBuf, fancy_ui: bool) -> Result<()> {
    use crate::model::{Alignment, FileBrowserState};

    let start_dir = if start_dir.is_dir() {
        start_dir
    } else {
        std::env::current_dir().unwrap_or_else(|_| std::path::PathBuf::from("."))
    };
    let mut state = AppState::new(Alignment::new(vec![]), "No file".to_string());
    state.fancy_ui = fancy_ui;
    state.file_browser = Some(FileBrowserState::new(start_dir, "Select a sequence file".to_string()));

    let mut app = App::new(state)?;
    app.run()
}

/// Convenience function to run the application with background loading.
/// Optional `preset_translation` is (genetic_code_id, reading_frame) to preset translation settings.
pub fn run_app_with_loading(
    file_path: PathBuf,
    forced_format: Option<FileFormat>,
    preset_translation: Option<(u8, u8)>,
    fancy_ui: bool,
) -> Result<()> {
    let mut app = App::new_with_background_load(file_path, forced_format, preset_translation, fancy_ui)?;
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
