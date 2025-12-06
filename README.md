# Rosalind Solutions for the Problem IDs 

All of the functions are basically implemented within a one file: `rosalind_functions.py`. 

## Implemented Problems

### 1. DNA – Counting DNA Nucleotides (`DNA`)
- **Function:** `count_dna_nucleotides(dna)`
- **Input:** DNA string (`A`, `C`, `G`, `T`)
- **Output:** Tuple `(A_count, C_count, G_count, T_count)`

### 2. RNA – Transcribing DNA into RNA (`RNA`)
- **Function:** `dna_to_rna(dna)`
- **Input:** DNA string
- **Output:** RNA string where `T` is replaced by `U`

### 3. REVC – Complementing a Strand of DNA (`REVC`)
- **Function:** `reverse_complement(dna)`
- **Input:** DNA string
- **Output:** Reverse complement DNA string

### 4. GC – Computing GC Content (`GC`)
- **Function:** `max_gc_content(fasta_str)`
- **Input:** Multi-FASTA formatted string
- **Output:** `(id, gc_percent)` of the sequence with the highest GC content

### 5. HAMM – Counting Point Mutations (`HAMM`)
- **Function:** `hamming_distance(s, t)`
- **Input:** Two equal-length DNA strings
- **Output:** Integer Hamming distance

### 6. PROT – Translating RNA into Protein (`PROT`)
- **Function:** `rna_to_protein(rna)`
- **Input:** RNA string
- **Output:** Protein string (translation stops at first stop codon)

### 7. SUBS – Finding a Motif in DNA (`SUBS`)
- **Function:** `motif_substring_locations(s, t)`
- **Input:** DNA string `s` and motif `t`
- **Output:** List of 1-based indices where `t` occurs in `s` (overlaps allowed)

### 8. PRTM – Calculating Protein Mass (`PRTM`)
- **Function:** `protein_mass(protein)`
- **Input:** Protein string
- **Output:** Total monoisotopic mass (float)

### 9. REVP – Locating Restriction Sites (`REVP`)
- **Function:** `reverse_palindromes(dna)`
- **Input:** DNA string
- **Output:** List of `(position, length)` for reverse palindromic sites (length 4–12)

### 10. TRAN – Transitions and Transversions (`TRAN`)
- **Function:** `transition_transversion_ratio(s1, s2)`
- **Input:** Two equal-length DNA strings
- **Output:** Ratio transitions / transversions (float)
