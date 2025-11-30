//! Genetic code definitions and translation logic.
//!
//! This module provides:
//! - NCBI genetic code tables (1-33)
//! - Codon to amino acid translation
//! - Support for different reading frames

/// A genetic code table for translating codons to amino acids.
#[derive(Debug, Clone)]
pub struct GeneticCode {
    /// NCBI genetic code ID
    pub id: u8,
    /// Name of the genetic code
    pub name: String,
    /// Codon to amino acid mapping as flat array (64 entries)
    /// Index = base1_idx * 16 + base2_idx * 4 + base3_idx
    /// where T=0, C=1, A=2, G=3
    codon_table: [u8; 64],
}

/// Convert nucleotide to index: T/U=0, C=1, A=2, G=3, invalid=255
#[inline(always)]
fn base_to_index(b: u8) -> u8 {
    match b {
        b'T' | b't' | b'U' | b'u' => 0,
        b'C' | b'c' => 1,
        b'A' | b'a' => 2,
        b'G' | b'g' => 3,
        _ => 255,
    }
}

/// Check if a nucleotide is a simple ambiguity code (R, Y, N, ?)
/// Returns the possible base indices, or empty slice if not ambiguous
#[inline(always)]
fn ambiguity_expansions(b: u8) -> &'static [u8] {
    match b {
        b'R' | b'r' => &[2, 3],     // A or G (puRine)
        b'Y' | b'y' => &[0, 1],     // T/U or C (pYrimidine)
        b'N' | b'n' | b'?' => &[0, 1, 2, 3], // any base
        _ => &[],
    }
}

impl GeneticCode {
    /// Creates a new genetic code from NCBI format strings.
    /// 
    /// # Arguments
    /// * `id` - NCBI genetic code ID
    /// * `name` - Name of the genetic code
    /// * `ncbieaa` - 64-character string of amino acids (NCBI format)
    fn new(id: u8, name: &str, ncbieaa: &str) -> Self {
        // NCBI order: TTT, TTC, TTA, TTG, TCT, TCC, ... (Base1, Base2, Base3)
        // This matches our index calculation: T=0, C=1, A=2, G=3
        let mut codon_table = [b'X'; 64];
        let aa_bytes = ncbieaa.as_bytes();
        for (i, &aa) in aa_bytes.iter().enumerate().take(64) {
            codon_table[i] = aa;
        }
        
        Self {
            id,
            name: name.to_string(),
            codon_table,
        }
    }
    
    /// Translates a single codon (3 bytes) to an amino acid.
    /// 
    /// # Rules:
    /// - Standard codons (3 nucleotides) are translated using the genetic code
    /// - Simple ambiguity codes (R, Y, N, ?) with 2 unambiguous bases:
    ///   If all possible codons translate to the same AA, return that AA
    /// - Codons with gaps (-) or frameshifts (!) follow special rules:
    ///   - All gaps/frameshifts (e.g., "---", "!!!") → '-'
    ///   - Mixed with nucleotides (e.g., "C-T", "A!G") → '!'
    /// - Other ambiguous/invalid codons → 'X'
    #[inline]
    pub fn translate_codon(&self, codon: &[u8]) -> u8 {
        if codon.len() != 3 {
            return b'X';
        }
        
        let b1 = codon[0];
        let b2 = codon[1];
        let b3 = codon[2];
        
        // Check for gaps and frameshifts
        let is_gap = |b: u8| b == b'-' || b == b'!' || b == b'.';
        let g1 = is_gap(b1);
        let g2 = is_gap(b2);
        let g3 = is_gap(b3);
        
        // All gaps/frameshifts → gap
        if g1 && g2 && g3 {
            return b'-';
        }
        
        // Mixed gaps/frameshifts with nucleotides → frameshift
        if g1 || g2 || g3 {
            return b'!';
        }
        
        // Convert to indices (handles case and U→T)
        let i1 = base_to_index(b1);
        let i2 = base_to_index(b2);
        let i3 = base_to_index(b3);
        
        // Fast path: all unambiguous bases
        if i1 != 255 && i2 != 255 && i3 != 255 {
            let idx = (i1 as usize) * 16 + (i2 as usize) * 4 + (i3 as usize);
            return self.codon_table[idx];
        }
        
        // Check for ambiguity codes (R, Y, N, ?)
        let a1 = ambiguity_expansions(b1);
        let a2 = ambiguity_expansions(b2);
        let a3 = ambiguity_expansions(b3);
        
        // Count ambiguous positions
        let amb_count = (i1 == 255 && !a1.is_empty()) as u8
            + (i2 == 255 && !a2.is_empty()) as u8
            + (i3 == 255 && !a3.is_empty()) as u8;
        
        // Only handle 1 ambiguous position (practical case)
        // Multiple ambiguities (e.g., NNR) → too many possibilities, return X
        if amb_count != 1 {
            return b'X';
        }
        
        // Get the expansion arrays (use single-element slices for unambiguous bases)
        let exp1: &[u8] = if i1 != 255 { std::slice::from_ref(&i1) } else { a1 };
        let exp2: &[u8] = if i2 != 255 { std::slice::from_ref(&i2) } else { a2 };
        let exp3: &[u8] = if i3 != 255 { std::slice::from_ref(&i3) } else { a3 };
        
        // Check if all possible codons translate to the same amino acid
        let mut first_aa: Option<u8> = None;
        for &e1 in exp1 {
            for &e2 in exp2 {
                for &e3 in exp3 {
                    let idx = (e1 as usize) * 16 + (e2 as usize) * 4 + (e3 as usize);
                    let aa = self.codon_table[idx];
                    match first_aa {
                        None => first_aa = Some(aa),
                        Some(prev) if prev != aa => return b'X', // Different AAs → X
                        _ => {}
                    }
                }
            }
        }
        
        first_aa.unwrap_or(b'X')
    }
    
    /// Translates a single codon string to an amino acid char (convenience wrapper).
    pub fn translate_codon_str(&self, codon: &str) -> char {
        self.translate_codon(codon.as_bytes()) as char
    }
    
    /// Translates an entire nucleotide sequence to amino acids.
    /// 
    /// # Arguments
    /// * `sequence` - The nucleotide sequence as bytes
    /// * `frame` - Reading frame (0, 1, or 2 for +1, +2, +3)
    /// 
    /// Returns the translated amino acids as a Vec<u8>.
    pub fn translate_sequence(&self, sequence: &[u8], frame: usize) -> Vec<u8> {
        let len = sequence.len();
        let start = frame.min(2);
        
        // Pre-allocate result with exact capacity
        let num_codons = (len.saturating_sub(start)) / 3;
        let mut result = Vec::with_capacity(num_codons);
        
        let mut pos = start;
        while pos + 3 <= len {
            result.push(self.translate_codon(&sequence[pos..pos+3]));
            pos += 3;
        }
        
        result
    }
    
    /// Translates a nucleotide sequence string to amino acids string.
    /// Convenience wrapper for translate_sequence.
    pub fn translate_sequence_str(&self, sequence: &str, frame: usize) -> String {
        let result = self.translate_sequence(sequence.as_bytes(), frame);
        // SAFETY: All amino acid chars are valid ASCII/UTF-8
        unsafe { String::from_utf8_unchecked(result) }
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
        assert_eq!(standard.translate_codon_str("ATG"), 'M'); // Start codon
        assert_eq!(standard.translate_codon_str("TAA"), '*'); // Stop codon
        assert_eq!(standard.translate_codon_str("TAG"), '*'); // Stop codon
        assert_eq!(standard.translate_codon_str("TGA"), '*'); // Stop codon
        assert_eq!(standard.translate_codon_str("TTT"), 'F'); // Phenylalanine
        assert_eq!(standard.translate_codon_str("GGG"), 'G'); // Glycine
    }
    
    #[test]
    fn test_rna_translation() {
        let codes = GeneticCodes::new();
        let standard = codes.default_code();
        
        // U should be treated as T
        assert_eq!(standard.translate_codon_str("AUG"), 'M');
        assert_eq!(standard.translate_codon_str("UUU"), 'F');
    }
    
    #[test]
    fn test_gap_handling() {
        let codes = GeneticCodes::new();
        let standard = codes.default_code();
        
        // All gaps → gap
        assert_eq!(standard.translate_codon_str("---"), '-');
        assert_eq!(standard.translate_codon_str("!!!"), '-');
        assert_eq!(standard.translate_codon_str("-!-"), '-');
        
        // Mixed → frameshift
        assert_eq!(standard.translate_codon_str("A--"), '!');
        assert_eq!(standard.translate_codon_str("-T-"), '!');
        assert_eq!(standard.translate_codon_str("AT-"), '!');
        assert_eq!(standard.translate_codon_str("A!G"), '!');
    }
    
    #[test]
    fn test_ambiguous_nucleotides() {
        let codes = GeneticCodes::new();
        let standard = codes.default_code();
        
        // R = A or G (purine)
        // CTR → CTA (Leu) or CTG (Leu) → both Leu, so L
        assert_eq!(standard.translate_codon_str("CTR"), 'L');
        
        // Y = T or C (pyrimidine)
        // TAY → TAT (Tyr) or TAC (Tyr) → both Tyr, so Y
        assert_eq!(standard.translate_codon_str("TAY"), 'Y');
        
        // N = any base
        // GGN → GGT, GGC, GGA, GGG → all Gly, so G
        assert_eq!(standard.translate_codon_str("GGN"), 'G');
        
        // ACN → ACT, ACC, ACA, ACG → all Thr, so T
        assert_eq!(standard.translate_codon_str("ACN"), 'T');
        
        // ? treated like N
        assert_eq!(standard.translate_codon_str("GG?"), 'G');
        
        // Cases where ambiguity leads to different AAs → X
        // ATN → ATT (Ile), ATC (Ile), ATA (Ile), ATG (Met) → mixed, X
        assert_eq!(standard.translate_codon_str("ATN"), 'X');
        
        // AGR → AGA (Arg), AGG (Arg) → both Arg, so R
        assert_eq!(standard.translate_codon_str("AGR"), 'R');
        
        // AGY → AGT (Ser), AGC (Ser) → both Ser, so S
        assert_eq!(standard.translate_codon_str("AGY"), 'S');
        
        // Multiple ambiguous positions → X (too complex)
        assert_eq!(standard.translate_codon_str("NNN"), 'X');
        assert_eq!(standard.translate_codon_str("NTN"), 'X');
        assert_eq!(standard.translate_codon_str("RRR"), 'X');
        
        // Unknown ambiguity codes → X
        assert_eq!(standard.translate_codon_str("ATW"), 'X'); // W = A or T
        assert_eq!(standard.translate_codon_str("ATS"), 'X'); // S = G or C
    }
    
    #[test]
    fn test_sequence_translation() {
        let codes = GeneticCodes::new();
        let standard = codes.default_code();
        
        let seq = "ATGTTTTAG"; // M, F, *
        assert_eq!(standard.translate_sequence_str(seq, 0), "MF*");
        
        // Different frames
        let seq = "AATGTTTTAG"; // 10 characters
        // Frame 0: AAT GTT TTA G → N, V, L (incomplete G not translated)
        assert_eq!(standard.translate_sequence_str(seq, 0), "NVL");
        // Frame 1: ATG TTT TAG → M, F, *
        assert_eq!(standard.translate_sequence_str(seq, 1), "MF*");
        // Frame 2: TGT TTT AG → C, F (incomplete AG not translated)
        assert_eq!(standard.translate_sequence_str(seq, 2), "CF");
    }
    
    #[test]
    fn test_case_insensitive() {
        let codes = GeneticCodes::new();
        let standard = codes.default_code();
        
        assert_eq!(standard.translate_codon_str("atg"), 'M');
        assert_eq!(standard.translate_codon_str("AtG"), 'M');
    }
    
    #[test]
    fn test_different_genetic_codes() {
        let codes = GeneticCodes::new();
        
        // In standard code, TGA is stop
        let standard = codes.get(1).unwrap();
        assert_eq!(standard.translate_codon_str("TGA"), '*');
        
        // In vertebrate mitochondrial (code 2), TGA is Trp (W)
        let vert_mito = codes.get(2).unwrap();
        assert_eq!(vert_mito.translate_codon_str("TGA"), 'W');
    }
}
