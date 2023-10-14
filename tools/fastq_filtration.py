EXAMPLE_FASTQ = {
    '@SRX079804:1:SRR292678:1:1101:21885:21885': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079804:1:SRR292678:1:1101:24563:24563': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079804:1:SRR292678:1:1101:30161:30161': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD'),
    '@SRX079804:1:SRR292678:1:1101:47176:47176': ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG', 'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB'),
    '@SRX079804:1:SRR292678:1:1101:149302:149302': ('TAGGGTTGTATTTGCAGATCCATGGCATGCCAAAAAGAACATCGTCCCGTCCAATATCTGCAACATACCAGTTGGTTGGTA', '@;CBA=:@;@DBDCDEEE/EEEEEEF@>FBEEB=EFA>EEBD=DAEEEEB9)99>B99BC)@,@<9CDD=C,5;B::?@;A'),
    '@SRX079804:1:SRR292678:1:1101:170868:170868': ('CTGCCGAGACTGTTCTCAGACATGGAAAGCTCGATTCGCATACACTCGCTGAGTAAGAGAGTCACACCAAATCACAGATT', 'E;FFFEGFGIGGFBG;C6D<@C7CDGFEFGFHDFEHHHBBHHFDFEFBAEEEEDE@A2=DA:??C3<BCA7@DCDEG*EB'),
    '@SRX079804:1:SRR292678:1:1101:171075:171075': ('CATTATAGTAATACGGAAGATGACTTGCTGTTATCATTACAGCTCCATCGCATGAATAATTCTCTAATATAGTTGTCAT', 'HGHHHHGFHHHHFHHEHHHHFGEHFGFGGGHHEEGHHEEHBHHFGDDECEGGGEFGF<FGGIIGEBGDFFFGFFGGFGF'),
    '@SRX079804:1:SRR292678:1:1101:175500:175500': ('GACGCCGTGGCTGCACTATTTGAGGCACCTGTCCTCGAAGGGAAGTTCATCTCGACGCGTGTCACTATGACATGAATG', 'GGGGGFFCFEEEFFDGFBGGGA5DG@5DDCBDDE=GFADDFF5BE49<<<BDD?CE<A<8:59;@C.C9CECBAC=DE'),
    '@SRX079804:1:SRR292678:1:1101:190136:190136': ('GAACCTTCTTTAATTTATCTAGAGCCCAAATTTTAGTCAATCTATCAACTAAAATACCTACTGCTACTACAAGTATT', 'DACD@BEECEDE.BEDDDDD,>:@>EEBEEHEFEHHFFHH?FGBGFBBD77B;;C?FFFFGGFED.BBABBG@DBBE'),
    '@SRX079804:1:SRR292678:1:1101:190845:190845': ('CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC', 'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE'),
    '@SRX079804:1:SRR292678:1:1101:198993:198993': ('AGTTATTTATGCATCATTCTCATGTATGAGCCAACAAGATAGTACAAGTTTTATTGCTATGAGTTCAGTACAACA', '<<<=;@B??@<>@><48876EADEG6B<A@*;398@.=BB<7:>.BB@.?+98204<:<>@?A=@EFEFFFEEFB'),
    '@SRX079804:1:SRR292678:1:1101:204480:204480': ('AGTGAGACACCCCTGAACATTCCTAGTAAGACATCTTTGAATATTACTAGTTAGCCACACTTTAAAATGACCCG', '<98;<@@@:@CD@BCCDD=DBBCEBBAAA@9???@BCDBCGF=GEGDFGDBEEEEEFFFF=EDEE=DCD@@BBC')
    } 

def filter_by_quality(seq_and_quality:tuple, quality_threshold:int) -> str:
    """
    Filter sequence by quality. Decoding quality sequence in q-score list and calculate quality. 

    Arguments:
    - seq_and_quality (tuple): pair of fastq-seq and quality sequence
    - quality_threshold (int): allow filter sequence and it is the bottom bound. Has got default value: 0. 

    Return:
    - str, if fastq-seq corresponds filtration condition.
    - None in opposite situation.
    """
    seq = seq_and_quality[0]
    quality_seq = seq_and_quality[1]

    q_score = []
    for symbol in quality_seq:
        q_score.append(ord(symbol)-33)

    quality = round(sum(q_score)/len(q_score), 4)
    if quality >= quality_threshold:
        return seq


def filtre_by_length(seq:str, length_bounds:tuple) -> str:
    """
    Filter sequence by length.

    Arguments:
    - seq (str): fastq-seq
    - length_bounds (tuple): bottom and top bound to filtration. Has got default value: (0, 2**32).

    Return:
    - str, if fastq-seq corresponds filtration condition.
    - None in opposite situation.
    """
    length_seq = len(seq)
    if length_seq >= length_bounds[0] and length_seq <= length_bounds[1]:
        return seq



def filter_by_gc_level(seq:str, gc_bounds:tuple) -> str:
    """
    Filter sequence by gc-level.

    Arguments:
    - seq (str): fastq-seq
    - gc_bounds (tuple): bottom and top bound to filtration. Use only top bound if gc_bound is int.

    Return:
    - str, if fastq-seq corresponds filtration condition.
    - None in opposite situation. 
    """
    if type(gc_bounds) == int:
        bottom_bound = 0
        top_bound = gc_bounds
    elif type(gc_bounds) == tuple:
        bottom_bound = gc_bounds[0]
        top_bound = gc_bounds[1]

    gc_amount = 0
    for nucl in seq:
        if nucl == "G" or nucl == "C":
            gc_amount += 1
    
    cur_gc_level = round(gc_amount / len(seq) * 100, 4)
    if cur_gc_level >= bottom_bound and cur_gc_level <= top_bound:
        return seq

def save_filtered_fastq(func_return:tuple) -> str:
    """
    
    """
    dict_of_filtered_fastq = func_return[0]
    output_filename = func_return[1]
    start_path = func_return[2]

    cur_dir = 'fastq_filtrator_resuls'

    if output_filename == "" and cur_dir not in os.listdir(os.getcwd()):
        output_filename = start_path
        os.mkdir(cur_dir)
    elif output_filename != "" and cur_dir not in os.listdir(os.getcwd()):
        os.mkdir(cur_dir)

    with open(os.path.join(cur_dir, output_filename), mode='w') as file:
        for fastq_seq in dict_of_filtered_fastq.items():
            file.write(fastq_seq[0]+'\n')
            file.write(fastq_seq[1][0]+'\n')
            file.write(fastq_seq[1][1]+'\n')
            file.write(fastq_seq[1][2]+'\n')
    return f"Fastq-seqs write in {output_filename} file."  

def filter_fastq(input_path:str, output_filename:str, gc_bounds:tuple, length_bounds=(0, 2**32), quality_threshold=0) -> dict:
    """
    Function containing conditions to filtration fastq-sequences.
    
    Takes arbitrary number of arguments with characteristic to filtrarion.
    Returns the result of the filtration in dictionary.

    Arguments:
    - dict_of_fastq (dict): dictionary with fastq-seq in format 'name' : ('sequence', 'quality')
    - gc_bounds (tuple or int): bottom and top bound to filtration. Use only top bound if gc_bound is int
    - length_bounds (tuple): bottom and top bound to filtration. Has got default value: (0, 2**32).
    - quality_threshold (int): allow filter sequence and it is the bottom bound. Has got default value: 0.
    """

    with open(os.path.abspath(input_path)) as fastq_file:

        fastq_lines = [line.strip('\n') for line in fastq_file.readlines()]

        dict_of_fastq = dict()
        for i in range(0, len(fastq_lines), 4):
            dict_of_fastq[fastq_lines[i]] = tuple(fastq_lines[i+1:i+4])


    seqs_with_qualities = [triplet for triplet in dict_of_fastq.values()]

    triplets = []
    for triplet in seqs_with_qualities:
        if (filter_by_gc_level(triplet[0], gc_bounds) and filtre_by_length(triplet[0], length_bounds) and filter_by_quality((triplet[0],triplet[2]), quality_threshold)) != None:
            triplets.append(triplet)

    ans_dict = {}
    for triplet in triplets:
        for item in dict_of_fastq.items():
            if triplet in item:
                ans_dict[item[0]] = triplet


    return ans_dict
