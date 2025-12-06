# DNA – Counting DNA Nucleotides
# Problem ID: DNA
def count_dna_nucleotides(dna: str):
    """
    Rosalind 'DNA' – Counting DNA Nucleotides

    Input: 
        dna (str): A DNA string consisting of characters 'A', 'C', 'G', 'T'.

    Return:
        tuple: (count_A, count_C, count_G, count_T)

    Example:
        >>> count_dna_nucleotides("AGCT")
        (1, 1, 1, 1)
    """
    a = c = g = t = 0
    for base in dna:
        if base == "A":
            a += 1
        elif base == "C":
            c += 1
        elif base == "G":
            g += 1
        elif base == "T":
            t += 1
    return a, c, g, t

    # RNA – Transcribing DNA into RNA
# Problem ID: RNA
def dna_to_rna(dna: str) -> str:
    """
    Rosalind 'RNA' – Transcribing DNA into RNA

    Input:
        dna (str): A DNA string with 'A', 'C', 'G', 'T'.

    Return:
        str: Corresponding RNA string (all 'T' replaced with 'U').
    """
    return dna.replace("T", "U")


# REVC – Complementing a Strand of DNA
# Problem ID: REVC
def reverse_complement(dna: str) -> str:
    """
    Rosalind 'REVC' – Complementing a Strand of DNA

    Input:
        dna (str): A DNA string.

    Return:
        str: Reverse complement of the input DNA.
    """
    comp = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }
    rev = []
    # reverse iterate
    for base in reversed(dna):
        rev.append(comp[base])
    return "".join(rev)


# GC – Computing GC Content
# Problem ID: GC
def max_gc_content(fasta_str: str):
    """
    Rosalind 'GC' – Computing GC Content

    Input:
        fasta_str (str): A multi-FASTA formatted string, e.g.:

            >Rosalind_6404
            CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTT
            >Rosalind_5959
            CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGAA...

    Return:
        tuple: (max_id, max_gc_percent)

            max_id (str): ID of sequence with highest GC content.
            max_gc_percent (float): GC percentage (0–100).

    Notes:
        - GC% is computed as 100 * (count(G)+count(C)) / length.
    """
    max_id = None
    max_gc = -1.0

    current_id = None
    current_seq = []

    lines = fasta_str.strip().splitlines()
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            # finalize previous sequence
            if current_id is not None:
                seq = "".join(current_seq)
                gc_count = seq.count("G") + seq.count("C")
                gc_percent = (gc_count / len(seq)) * 100.0 if seq else 0.0
                if gc_percent > max_gc:
                    max_gc = gc_percent
                    max_id = current_id
            # start new sequence
            current_id = line[1:]  # drop '>'
            current_seq = []
        else:
            current_seq.append(line)

    # finalize last seq
    if current_id is not None:
        seq = "".join(current_seq)
        gc_count = seq.count("G") + seq.count("C")
        gc_percent = (gc_count / len(seq)) * 100.0 if seq else 0.0
        if gc_percent > max_gc:
            max_gc = gc_percent
            max_id = current_id

    return max_id, max_gc


# HAMM – Counting Point Mutations
# Problem ID: HAMM
def hamming_distance(s: str, t: str) -> int:
    """
    Rosalind 'HAMM' – Counting Point Mutations

    Input:
        s (str), t (str): Two DNA strings of equal length.

    Return:
        int: The Hamming distance (number of differing positions).
    """
    # assumes len(s) == len(t)
    distance = 0
    for c1, c2 in zip(s, t):
        if c1 != c2:
            distance += 1
    return distance


# PROT – Translating RNA into Protein
# Problem ID: PROT
def rna_to_protein(rna: str) -> str:
    """
    Rosalind 'PROT' – Translating RNA into Protein

    Input:
        rna (str): An RNA string.

    Return:
        str: The translated protein sequence (stops at first stop codon).
    """
    codon_table = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "UAU": "Y", "UAC": "Y", "UAA": "Stop", "UAG": "Stop",
        "UGU": "C", "UGC": "C", "UGA": "Stop", "UGG": "W",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    protein = []
    # read in steps of 3
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        aa = codon_table.get(codon, "")
        if aa == "Stop":
            break
        if aa:  # ignore unknown codons
            protein.append(aa)
    return "".join(protein)


# SUBS – Finding a Motif in DNA
# Problem ID: SUBS
def motif_substring_locations(s: str, t: str):
    """
    Rosalind 'SUBS' – Finding a Motif in DNA

    Input:
        s (str): DNA string (text).
        t (str): DNA string (motif).

    Return:
        list[int]: 1-based starting positions of all (possibly overlapping)
                   occurrences of t in s.

    Example:
        >>> motif_substring_locations("GATATATGCATATACTT", "ATAT")
        [2, 4, 10]
    """
    positions = []
    len_s = len(s)
    len_t = len(t)
    for i in range(len_s - len_t + 1):
        if s[i:i+len_t] == t:
            positions.append(i + 1)  # 1-based
    return positions
