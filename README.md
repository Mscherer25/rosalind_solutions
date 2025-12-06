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
