//! Genetic code definitions and translation logic.
//!
//! This module provides:
//! - NCBI genetic code tables (1-33)
//! - Codon to amino acid translation
//! - Support for different reading frames

use std::collections::HashMap;

/// A genetic code table for translating codons to amino acids.
#[derive(Debug, Clone)]
pub struct GeneticCode {
    /// NCBI genetic code ID
    pub id: u8,
    /// Name of the genetic code
    pub name: String,
    /// Codon to amino acid mapping (64 entries for standard codons)
    codon_table: HashMap<String, char>,
}

impl GeneticCode {
    /// Creates a new genetic code from NCBI format strings.
    /// 
    /// # Arguments
    /// * `id` - NCBI genetic code ID
    /// * `name` - Name of the genetic code
    /// * `ncbieaa` - 64-character string of amino acids (NCBI format)
    fn new(id: u8, name: &str, ncbieaa: &str) -> Self {
        let bases = ['T', 'C', 'A', 'G'];
        let mut codon_table = HashMap::new();
        
        // NCBI order: TTT, TTC, TTA, TTG, TCT, TCC, ... (Base1, Base2, Base3)
        let mut idx = 0;
        for &b1 in &bases {
            for &b2 in &bases {
                for &b3 in &bases {
                    let codon = format!("{}{}{}", b1, b2, b3);
                    let aa = ncbieaa.chars().nth(idx).unwrap_or('X');
                    codon_table.insert(codon, aa);
                    idx += 1;
                }
            }
        }
        
        Self {
            id,
            name: name.to_string(),
            codon_table,
        }
    }
    
    /// Translates a single codon to an amino acid.
    /// 
    /// # Rules:
    /// - Standard codons (3 nucleotides) are translated using the genetic code
    /// - Codons with ambiguous nucleotides (not A, C, G, T/U) return 'X'
    /// - Codons with gaps (-) or frameshifts (!) follow special rules:
    ///   - All gaps/frameshifts (e.g., "---", "!!!") → '-'
    ///   - Mixed with nucleotides (e.g., "C-T", "A!G") → '!'
    pub fn translate_codon(&self, codon: &str) -> char {
        if codon.len() != 3 {
            return 'X';
        }
        
        let codon_upper: String = codon.to_uppercase();
        let chars: Vec<char> = codon_upper.chars().collect();
        
        // Check for gaps and frameshifts
        let gap_or_frameshift: Vec<bool> = chars.iter()
            .map(|&c| c == '-' || c == '!' || c == '.')
            .collect();
        
        // All gaps/frameshifts → gap
        if gap_or_frameshift.iter().all(|&b| b) {
            return '-';
        }
        
        // Mixed gaps/frameshifts with nucleotides → frameshift
        if gap_or_frameshift.iter().any(|&b| b) {
            return '!';
        }
        
        // Convert U to T for RNA
        let codon_dna: String = chars.iter()
            .map(|&c| if c == 'U' { 'T' } else { c })
            .collect();
        
        // Check for valid nucleotides
        let valid_nucleotides = codon_dna.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T'));
        
        if !valid_nucleotides {
            return 'X'; // Ambiguous nucleotide
        }
        
        // Look up in codon table
        self.codon_table.get(&codon_dna).copied().unwrap_or('X')
    }
    
    /// Translates an entire nucleotide sequence to amino acids.
    /// 
    /// # Arguments
    /// * `sequence` - The nucleotide sequence to translate
    /// * `frame` - Reading frame (0, 1, or 2 for +1, +2, +3)
    pub fn translate_sequence(&self, sequence: &str, frame: usize) -> String {
        let chars: Vec<char> = sequence.chars().collect();
        let mut result = String::new();
        
        let start = frame.min(2);
        let mut pos = start;
        
        while pos + 3 <= chars.len() {
            let codon: String = chars[pos..pos+3].iter().collect();
            result.push(self.translate_codon(&codon));
            pos += 3;
        }
        
        result
    }
}

/// All available genetic codes from NCBI.
pub struct GeneticCodes {
    codes: Vec<GeneticCode>,
}

impl GeneticCodes {
    /// Creates the complete set of NCBI genetic codes.
    pub fn new() -> Self {
        let codes = vec![
            GeneticCode::new(1, "Standard", 
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(2, "Vertebrate Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"),
            GeneticCode::new(3, "Yeast Mitochondrial",
                "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(4, "Mold/Protozoan/Coelenterate Mito...",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(5, "Invertebrate Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"),
            GeneticCode::new(6, "Ciliate/Dasycladacean/Hexamita Nuclear",
                "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(9, "Echinoderm/Flatworm Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
            GeneticCode::new(10, "Euplotid Nuclear",
                "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(11, "Bacterial/Archaeal/Plant Plastid",
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(12, "Alternative Yeast Nuclear",
                "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(13, "Ascidian Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"),
            GeneticCode::new(14, "Alternative Flatworm Mitochondrial",
                "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
            GeneticCode::new(15, "Blepharisma Macronuclear",
                "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(16, "Chlorophycean Mitochondrial",
                "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(21, "Trematode Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
            GeneticCode::new(22, "Scenedesmus obliquus Mitochondrial",
                "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(23, "Thraustochytrium Mitochondrial",
                "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(24, "Rhabdopleuridae Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"),
            GeneticCode::new(25, "Candidate Division SR1/Gracilibacteria",
                "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(26, "Pachysolen tannophilus Nuclear",
                "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(27, "Karyorelict Nuclear",
                "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(28, "Condylostoma Nuclear",
                "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(29, "Mesodinium Nuclear",
                "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(30, "Peritrich Nuclear",
                "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(31, "Blastocrithidia Nuclear",
                "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(32, "Balanophoraceae Plastid",
                "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
            GeneticCode::new(33, "Cephalodiscidae Mitochondrial",
                "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"),
        ];
        
        Self { codes }
    }
    
    /// Returns all genetic codes.
    pub fn all(&self) -> &[GeneticCode] {
        &self.codes
    }
    
    /// Gets a genetic code by ID.
    pub fn get(&self, id: u8) -> Option<&GeneticCode> {
        self.codes.iter().find(|c| c.id == id)
    }
    
    /// Gets the default (Standard) genetic code.
    pub fn default_code(&self) -> &GeneticCode {
        self.get(1).expect("Standard genetic code should always exist")
    }
}

impl Default for GeneticCodes {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_standard_code_translation() {
        let codes = GeneticCodes::new();
        let standard = codes.default_code();
        
        // Test some common codons
        assert_eq!(standard.translate_codon("ATG"), 'M'); // Start codon
        assert_eq!(standard.translate_codon("TAA"), '*'); // Stop codon
        assert_eq!(standard.translate_codon("TAG"), '*'); // Stop codon
        assert_eq!(standard.translate_codon("TGA"), '*'); // Stop codon
        assert_eq!(standard.translate_codon("TTT"), 'F'); // Phenylalanine
        assert_eq!(standard.translate_codon("GGG"), 'G'); // Glycine
    }
    
    #[test]
    fn test_rna_translation() {
        let codes = GeneticCodes::new();
        let standard = codes.default_code();
        
        // U should be treated as T
        assert_eq!(standard.translate_codon("AUG"), 'M');
        assert_eq!(standard.translate_codon("UUU"), 'F');
    }
    
    #[test]
    fn test_gap_handling() {
        let codes = GeneticCodes::new();
        let standard = codes.default_code();
        
        // All gaps → gap
        assert_eq!(standard.translate_codon("---"), '-');
        assert_eq!(standard.translate_codon("!!!"), '-');
        assert_eq!(standard.translate_codon("-!-"), '-');
        
        // Mixed → frameshift
        assert_eq!(standard.translate_codon("A--"), '!');
        assert_eq!(standard.translate_codon("-T-"), '!');
        assert_eq!(standard.translate_codon("AT-"), '!');
        assert_eq!(standard.translate_codon("A!G"), '!');
    }
    
    #[test]
    fn test_ambiguous_nucleotides() {
        let codes = GeneticCodes::new();
        let standard = codes.default_code();
        
        // Ambiguous nucleotides → X
        assert_eq!(standard.translate_codon("ATN"), 'X');
        assert_eq!(standard.translate_codon("NNN"), 'X');
        assert_eq!(standard.translate_codon("CTR"), 'X'); // R = A or G
    }
    
    #[test]
    fn test_sequence_translation() {
        let codes = GeneticCodes::new();
        let standard = codes.default_code();
        
        let seq = "ATGTTTTAG"; // M, F, *
        assert_eq!(standard.translate_sequence(seq, 0), "MF*");
        
        // Different frames
        let seq = "AATGTTTTAG"; // 10 characters
        // Frame 0: AAT GTT TTA G → N, V, L (incomplete G not translated)
        assert_eq!(standard.translate_sequence(seq, 0), "NVL");
        // Frame 1: ATG TTT TAG → M, F, *
        assert_eq!(standard.translate_sequence(seq, 1), "MF*");
        // Frame 2: TGT TTT AG → C, F (incomplete AG not translated)
        assert_eq!(standard.translate_sequence(seq, 2), "CF");
    }
    
    #[test]
    fn test_case_insensitive() {
        let codes = GeneticCodes::new();
        let standard = codes.default_code();
        
        assert_eq!(standard.translate_codon("atg"), 'M');
        assert_eq!(standard.translate_codon("AtG"), 'M');
    }
    
    #[test]
    fn test_different_genetic_codes() {
        let codes = GeneticCodes::new();
        
        // In standard code, TGA is stop
        let standard = codes.get(1).unwrap();
        assert_eq!(standard.translate_codon("TGA"), '*');
        
        // In vertebrate mitochondrial (code 2), TGA is Trp (W)
        let vert_mito = codes.get(2).unwrap();
        assert_eq!(vert_mito.translate_codon("TGA"), 'W');
    }
}
