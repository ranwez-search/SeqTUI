// Use jemalloc
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use std::time::Instant;
use std::fs;

fn main() {
    let path = "/Users/ranwez/Downloads/triticeae_allindividuals_OneCopyGenes.fasta";
    
    println!("=== PROFILING WITH JEMALLOC ===");
    print_memory("Before loading");
    
    // Load file
    let start = Instant::now();
    let content = fs::read_to_string(path).unwrap();
    println!("File read: {:?}, size: {} MB", start.elapsed(), content.len() / 1_000_000);
    print_memory("After file read");
    
    // Parse
    let start = Instant::now();
    let alignment = seqtui::formats::fasta::parse_fasta_str(&content).unwrap();
    println!("Parse: {:?}, {} sequences, {} length", 
             start.elapsed(), 
             alignment.sequence_count(),
             alignment.alignment_length());
    drop(content);
    print_memory("After parse (file dropped)");
    
    // Force loading all sequence data
    println!("\n--- Force loading all NT data ---");
    let start = Instant::now();
    let mut total_nt = 0u64;
    for seq in &alignment.sequences {
        for &b in seq.as_bytes() {
            total_nt += b as u64;
        }
    }
    println!("Read all NT: {:?}, checksum: {}", start.elapsed(), total_nt);
    print_memory("After reading all NT");
    
    // Translation
    let codes = seqtui::genetic_code::GeneticCodes::new();
    let code = codes.default_code();
    
    println!("\n--- Translation ---");
    let start = Instant::now();
    let mut total_aa = 0usize;
    let translated: Vec<Vec<u8>> = alignment.sequences.iter()
        .map(|seq| {
            let aa = code.translate_sequence(seq.as_bytes(), 0);
            total_aa += aa.len();
            aa
        })
        .collect();
    println!("Translation: {:?}, total AA: {} MB", start.elapsed(), total_aa / 1_000_000);
    print_memory("After translation");
    
    println!("\nSequence data size: {} MB", 
        alignment.sequences.iter().map(|s| s.len()).sum::<usize>() / 1_000_000);
    println!("Translated data size: {} MB", 
        translated.iter().map(|s| s.len()).sum::<usize>() / 1_000_000);
    
    // Drop translated and see if memory goes down
    println!("\n--- Drop translated ---");
    drop(translated);
    print_memory("After dropping translated");
    
    std::thread::sleep(std::time::Duration::from_secs(2));
    print_memory("Final (after 2s)");
}

fn print_memory(label: &str) {
    let pid = std::process::id();
    let output = std::process::Command::new("ps")
        .args(["-o", "rss=", "-p", &pid.to_string()])
        .output()
        .ok();
    
    if let Some(out) = output {
        let rss = String::from_utf8_lossy(&out.stdout);
        let rss_kb: u64 = rss.trim().parse().unwrap_or(0);
        println!("  {} --> RSS = {} MB", label, rss_kb / 1024);
    }
}
