import os
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
    Save filtered fastq-sequencies and their characteristics in directory,
    which called fastq_filtrator_resuls.
    If directory is absent, create it. If output_filename is empty, use start_path to name file.

    Arguments:
    - dict_of_filtered_fastq (dict): dict with filtered fastq-seqs
    - start_path (str): name of starting file
    - output_filename (str): name of output file
    
    Return None
    """
    if output_filename != "" and output_filename.split(".")[1] != "fastq":
        raise ValueError(f"Not a FASTQ format of {output_filename} file")

    cur_dir = 'fastq_filtrator_resuls'

    if output_filename == "" and cur_dir not in os.listdir(os.getcwd()):
        output_filename = start_path
        os.mkdir(cur_dir)
    elif output_filename != "" and cur_dir not in os.listdir(os.getcwd()):
        os.mkdir(cur_dir)
    elif output_filename == "" and cur_dir in os.listdir(os.getcwd()):
        output_filename = start_path

    with open(os.path.join(cur_dir, output_filename), mode='w') as file:
        for fastq_seq in dict_of_filtered_fastq.items():
            file.write(fastq_seq[0]+'\n')
            file.write(fastq_seq[1][0]+'\n')
            file.write(fastq_seq[1][1]+'\n')
            file.write(fastq_seq[1][2]+'\n') 

def read_fastq_file(input_path:str) -> dict:
    """
    Read lines in fastq-file.
    
    Arguments:
    - input_path (str): name of file with fastq-seqs

    Return dictionary in the format 'name':('seq', '')
    """
    with open(os.path.abspath(input_path)) as fastq_file:

        fastq_lines = [line.strip('\n') for line in fastq_file.readlines()]

        dict_of_fastq = dict()
        for i in range(0, len(fastq_lines), 4):
            dict_of_fastq[fastq_lines[i]] = tuple(fastq_lines[i+1:i+4])
    
    return dict_of_fastq

def filter_fastq(input_path:str, gc_bounds:tuple, length_bounds=(0, 2**32), quality_threshold=0, output_filename="") -> str:
    """
    Function containing conditions to filtration fastq-sequences.
    
    Takes arbitrary number of arguments with characteristic to filtrarion.
    Returns the result of the filtration in dictionary.

    Arguments:
    - input_path (str): name of file with fastq-seqs
    - dict_of_fastq (dict): dictionary with fastq-seq in format {'name' : ('sequence', 'comment', 'quality)}
    - gc_bounds (tuple or int): bottom and top bound to filtration. Use only top bound if gc_bound is int
    - length_bounds (tuple): bottom and top bound to filtration. Has got default value: (0, 2**32).
    - quality_threshold (int): allow filter sequence and it is the bottom bound. Has got default value: 0.
    - output_filename (str): file, where writes filtered fastq-seqs
    """
    dict_of_fastq = read_fastq_file(input_path)
    
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

    save_filtered_fastq(ans_dict, input_path, output_filename)
    
    if output_filename == "":
        return  f"FASTQ-seq write in {input_path} file."
    return f"FASTQ-seq write in {output_filename} file."
