DNA_ALPHABET = {"A", "T", "G", "C", "a", "t", "g", "c"}
RNA_ALPHABET = {"A", "U", "G", "C", "a", "u", "g", "c"}
RNA_COMPLEMENT_NUCLS = { "A":"U", "U":"A", "a":"u", "u":"a", "G":"C", "C":"G", "c":"g", "g":"c"}
DNA_COMPLEMENT_NUCLS = {"A":"T", "T":"A", "a": "t", "t":"a", "G":"C", "C":"G", "g":"c", "c":"g"}

def check_nucl_acid(seq:str, nucl_type) -> bool:
    """
    Check the sequence is RNA or DNA, return boolean.
    """
    if nucl_type == "rna":
        possible_nucls = RNA_ALPHABET
    elif nucl_type == "dna":
        possible_nucls = DNA_ALPHABET

    unique_chars = set(seq.upper())
    return unique_chars <= possible_nucls
    
def reverse_complement(seq:str) -> str:
    """
    Perform reverse complementary complementary sequence of RNA or DNA.
    
    Arguments:
    - seq (str): DNA or RNA sequence in upper or lower case.

    Return:
    - If input DNA or RNA return reverse complementary sequence.
    - If input rubbish return message "It isn't RNA or DNA."
    
    Examples:
    - "ATCATAGA","CAGatAGC" -> ['TCTATGAT', 'GCTatCTG']
    - "ACACARRR","UATGC","UUuAUAUA" -> ["It isn't RNA or DNA.", "It isn't RNA or DNA.", 'UAUAUaAA']
    """
    if check_nucl_acid(seq, "dna") or check_nucl_acid(seq, "rna"):
        return reverse(complement(seq))    
    raise ValueError(f"{seq} isn't RNA or DNA.")    
    
def complement(seq:str) -> str:
    """
    Perform complementary complementary sequence of RNA or DNA.

    Arguments:
    - seq (str): DNA or RNA sequence in upper or lower case.

    Return:
    - If input DNA or RNA return complementary sequence.
    - If input rubbish return message "It isn't RNA or DNA."

    Examples:
    - "ATCATAGA","CAGatAGC" -> ['TAGTATCT', 'GTCtaTCG']
    - "ACACARRR","UATGC","UUuAUAUA" -> ["It isn't RNA or DNA.", "It isn't RNA or DNA.", 'AAaUAUAU'
    """    
    compl_nucls_map = {}
    if check_nucl_acid(seq, "dna"):
        compl_nucls_map = DNA_COMPLEMENT_NUCLS
    elif check_nucl_acid(seq, "rna"):
        compl_nucls_map = RNA_COMPLEMENT_NUCLS
    else:
        raise ValueError(f"{seq} isn't RNA or DNA.")
    
    trans_seq = ""
    for nucl in seq:
        trans_seq += compl_nucls_map[nucl]
    return trans_seq
    
def reverse(seq:str) -> str:
    """
    Perform reverse DNA or RNA sequence.

    Arguments:
    - seq (str): DNA or RNA sequence in upper or lower case.

    Return:
    - If input DNA or RNA return reverse sequence.
    - If input rubbish return message "It isn't RNA or DNA."

    Examples:
    - "ATCATAGA","CAGatAGC" -> ['AGATACTA', 'CGAtaGAC']
    - "ACACARRR","UATGC","UUUAUAUA" -> ["It isn't RNA or DNA.", "It isn't RNA or DNA.", 'AUAUAUUU']
    """
    if check_nucl_acid(seq, "dna") or check_nucl_acid(seq, "rna"):
        return seq[::-1]
    raise ValueError(f"{seq} isn't RNA or DNA.")
       
def transcribe_dna(seq:str) -> str:
    """
    Perform the transcribtion of DNA into mRNA.

     Argument:
    - seq (str): DNA sequence in upper or lower case. 

    Return:
    - If input DNA sequence return mRNA.
    - If input RNA sequence return message "seq is RNA."
    - If input rubbosh return message "It isn't RNA or DNA."

    Examples:
    - "ATCATAGA","CAGatAGC" -> ['AUCAUAGA', 'CAGauAGC']
    - "ACACARRR","UATGC","UUUAUAUA" -> ["It isn't RNA or DNA", "It isn't RNA or DNA", 'UUUAUAUA is RNA.']
    """
    if check_nucl_acid(seq, "dna"):
        return seq.replace("T", "U").replace("t", "u")
    elif check_nucl_acid(seq, "rna"):
        raise ValueError(f"{seq} is RNA.")
    raise ValueError(f"{seq} isn't RNA or DNA.")
    
def run_dna_rna_tools(*args) -> list: 
    """
    Function containing methods for nucleic acids analysis.

    Takes arbitrary number of arguments with nucleic acids sequencies
    and the name of the procedure to be performed (always the last
    argument). Returns the result of the procedure as string or
    list of strings if input more than one sequence.
    """
    *seqs, procedure = args
    
    operations = {
        "transcribe_dna":transcribe_dna, 
        "reverse":reverse, 
        "complement":complement, 
        "reverse_complement":reverse_complement
    }
    new_seqs = []
    for seq in seqs:
        new_seqs.append(operations[procedure](seq))
    
    if len(new_seqs) == 1:
        return new_seqs[0] 
    return new_seqs
