ALPHABET_PROTEIN = {"A","C","D","E","F","G","H","I","K","L","M",
                    "N","P","Q","R","S","T","V","W","Y",}
ALPHABET_RNA = {"A", "U", "G", "C"}
AMINO_ACID_MASS = {
    "A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694, "C": 103.00919, "Q": 128.05858,"E": 129.04259,
    "G": 57.02146, "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,"M": 131.04049, "F": 147.06841,
    "P": 97.05276, "S": 87.03203, "T": 101.04768,"W": 186.07931,"Y": 163.06333,"V": 99.06841,
}
GYDROPHOBIC_AMINOACIDS = {"A", "V", "L", "I", "P", "F", "W", "M"}
DNA_CODONS = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "F": ["TTT", "TTC"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"],
    "K": ["AAA", "AAG"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "M": ["ATG"],
    "N": ["AAT", "AAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "W": ["TGG"],
    "Y": ["TAT", "TAC"],
    "*": ["UAA", "UAG", "UGA"],
}
RNA_CODONS = {
    "F": ["UUC", "UUU"],
    "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
    "I": ["AUU", "AUC", "AUA"],
    "M": ["AUG"],
    "V": ["GUU", "GUC", "GUA", "GUG"],
    "S": ["UCU", "UCC", "UCA", "UCG"],
    "P": ["CCU", "CCC", "CCA", "CCG"],
    "T": ["ACU", "ACC", "ACA", "ACG"],
    "A": ["GCU", "GCC", "GCA", "GCG"],
    "Y": ["UAC", "UAU"],
    "*": ["UAA", "UAG", "UGA"],
    "H": ["CAU", "CAC"],
    "Q": ["CAA", "CAG"],
    "N": ["AAU", "AAC"],
    "K": ["AAA", "AAG"],
    "D": ["GAU", "GAC"],
    "E": ["GAA", "GAG"],
    "C": ["UGU", "UGC"],
    "W": ["UGG"],
    "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "S": ["AGU", "AGC"],
    "G": ["GGU", "GGC", "GGA", "GGG"],
}

def is_protein(seq: str):
    """
    Check the sequence is protein return boolean.
    """
    unique_chars = set(seq.upper())
    return unique_chars <= ALPHABET_PROTEIN

def is_rna(seq: str):
    """
    Check the sequence is RNA, return boolean.
    """
    unique_chars = set(seq.upper())
    return unique_chars <= ALPHABET_RNA

def compute_molecular_weight(protein: str) -> tuple:
    """
    Compute molecular weight (g/mol) of protein sequence.

    Argument:
    - protein (str): protein sequence.

    Return:
    - tuple with protein sequence and computed molecular
    weight (float rounded to 3 decimal places).
    """
    molecular_weight = 0
    for amino_acid in protein.upper():
        molecular_weight += AMINO_ACID_MASS[amino_acid]
    return protein, round(molecular_weight, 3)

def compute_length(protein: str) -> tuple:
    """
    Compute the length of the input protein sequence.

     Argument:
    - protein (str): protein sequence.

    Return:
    - tuple with protein sequence and computed length.
    """
    return protein, len(protein)

def protein_to_dna(protein: str) -> str:
    """
    Returns possible variants of DNAs for a given protein sequence.

    Argument:
    - protein (str): protein sequence.

    Return:
    - string, variants of nucleic acids.
    If several codons correspond to a given amino acid they are displayed with a '/'.

    Does not distinguish between lowercase and uppercase letters.

    Examples:

    -'MACDRS' -> 'ATG GCT/GCC/GCA/GCG TGT/TGC GAT/GAC CGT/CGC/CGA/CGG/AGA/AGG TCT/TCC/TCA/TCG/AGT/AGC'
    -'MaCdrS' -> 'ATG GCT/GCC/GCA/GCG TGT/TGC GAT/GAC CGT/CGC/CGA/CGG/AGA/AGG TCT/TCC/TCA/TCG/AGT/AGC'
    """
    nucleic_acid_seq = ""

    for aa in protein.upper():
        codons = DNA_CODONS[aa]
        nucleic_acid_seq += "/".join(codons) + " "

    return nucleic_acid_seq[:-1]

def count_amino_acids(protein: str) -> dict:
    """
    Calculates the number of each aminoacid in a given protein sequence.

    Argument:
    - protein (str): protein sequence.

    Return:
    - dictionary, where a key is the aminoacid letter and value is number of this aminoacid.

    Does not distinguish between lowercase and uppercase letters.

    Examples:

    -'MACDRS' -> {'M': 1, 'A': 1, 'C': 1, 'D': 1, 'R': 1, 'S': 1}
    -'MaCdrS' -> {'M': 1, 'A': 1, 'C': 1, 'D': 1, 'R': 1, 'S': 1}
    """
    amino_acids_dict = {}
    for aa in protein.upper():
        if aa in amino_acids_dict:
            amino_acids_dict[aa] += 1
        else:
            amino_acids_dict[aa] = 1
    return amino_acids_dict

def compute_hydrophobicity(protein: str) -> tuple:
    """
    Compute the percentage of gydrophobic aminoacids in protein sequence.

    Argument:
    - protein (str): protein sequence. Includes hydrophobic
    and hydrophilic aminoacids.

    Return:
    - tuple with protein sequence and computed percentage
    of gydrophobic aminoacids.
    """
    count_of_gydrophobic = 0
    for aa in protein:
        if aa in GYDROPHOBIC_AMINOACIDS:
            count_of_gydrophobic += 1

    percentage = round(count_of_gydrophobic / len(protein) * 100, 3)

    return protein, percentage

def translate_rna(rna: str) -> str:
    """
    Perform the translation of mRNA seguence into protein sequence.

    Argument:
    - rna (str): mRNA sequence. Must contain start-codon and one of
    the stop-codons.

    Return:
    - str, protein sequence after translation.
    Always starts with "M" and ends with "*".
    """
    triplets = [rna[i : i + 3].upper() for i in range(0, len(rna), 3)]
    protein = []
    for triplet in triplets:
        for aminoacid in RNA_CODONS.keys():
            if triplet in RNA_CODONS[aminoacid]:
                protein.append(aminoacid)

    if protein[-1] != "*":
        raise ValueError("Stop-codon (*) is absent in mRNA")
    if protein[0] != "M":
        raise ValueError("Start-codon (M) is absent in mRNA")

    start = protein.index("M")
    stop = protein.index("*")
    return "".join(protein[start : stop + 1])

def check_mutations(rna: str, protein: str) -> str:
    """
    Check missense mutations in the protein sequence after translation.

    Uses additional function "translate_rna(seq)".

    Arguments:
    - rna (str): sequence of mRNA with/without mutations.
    Must contain start-codon and one of the stop-codons.
    - protein (str): protein sequence translated from mRNA.
    Must start with "M" and ends with "*" (stop-codon).

    Note: is_protein(seq) doesn't see "*", but it's used in the other part of function.

    Return:
    - str, if mRNA without mutations return "Protein without mutations."
    If there are mutations in protein, returns aminoacid(s) and their position(s)

    Examples:
    - "AUGGUAGGGAAAUUUUGA", "MVGKF*" ->  "Protein without mutations."
    - "AUGGUAGGGAAAUUUUGA", "MGGVF*" -> "Mutations:G2, V4."
    - "AUGGUAGGGAAAUUUUGA", "MGGKF" –> "ValueError: Stop (*) is absent"
    - "AUGGUAGGGAAAUUUUGA", "GGKF*" –> "ValueError: Start (M) is absent"
    - "AUGAAAAAAUGA", "MK*" -> "ValueError: Different length of translated protein and protein"
    """
    correct_protein = translate_rna(rna)
    bank_of_mutations = []

    if is_protein(protein[:-1]) is not True:
        raise ValueError("Invalid protein sequence")
    if is_rna(rna) is not True:
        raise ValueError("Invalid RNA sequence")
    if protein[-1] != "*":
        raise ValueError("Stop (*) is absent")
    if protein[0] != "M":
        raise ValueError("Start (M) is absent")
    if len(protein) != len(rna) / 3:
        raise ValueError("Different length of translated protein and protein")

    for i in range(len(correct_protein)):
        if correct_protein[i] != protein[i]:
            bank_of_mutations.append(f"{protein[i]}{i + 1}")

    if len(bank_of_mutations) == 0:
        return "Protein without mutations."
    else:
        return "Mutations: " + ", ".join(bank_of_mutations) + "."

def run_protein_tools(*args: str):
    """
    Function containing methods for protein analysis.

    Takes arbitrary number of arguments with protein sequencies
    and the name of the procedure to be performed (always the last
    argument). Returns the result of the procedure as string, tuple
    or dictionary if one sequnce is submitted or list if several.

    Note: if procedure 'check_mutations' is used then input must
    contain only three arguments: RNA sequence, protein sequence
    and the name of procedure itself.
    """
    *seqs, operation = args
    results = []
    operations = {
        "compute_molecular_weight": compute_molecular_weight,
        "compute_length": compute_length,
        "compute_hydrophobicity": compute_hydrophobicity,
        "count_amino_acids": count_amino_acids,
        "protein_to_dna": protein_to_dna

    }
    if operation == "check_mutations":
        results.append(check_mutations(seqs[0], seqs[1]))
    else:
        for seq in seqs:
            if is_protein(seq) is not True:
                raise ValueError(f"Invalid protein sequence ({seq})")
            if operation not in operations:
                raise ValueError(f"Wrong operation name ({operation})")
            else:
                results.append(operations[operation](seq))
    if len(results) == 1:
        return results[0]
    else:
        return results
