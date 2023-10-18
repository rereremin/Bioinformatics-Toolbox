import os

def convert_multiline_fasta_to_oneline(input_fasta:str, output_fasta="") -> str:
    """
    Converts multiline fasta to oneline and writes in fasta file. 
    If function works correct, it prints message "FASTA-seq write in {output_file} file".
    If output_fasts is empty, raise ValueError("File hasn't got name!").

    Arguments:
    - input_fasta (str): multiline fasta file 
    - output_fasta (str): name of file, where saved lines, has default value equal ""
    """
    with open(os.path.abspath(input_fasta)) as fasta_file:
        
        list_of_lines = []
        for line in fasta_file.readlines():
            if line.startswith(">"):
                line = "\t"
            list_of_lines.append(line)

    longlines_to_file = "".join(list_of_lines[1:]).replace("\n", "").replace("\t", "\n")

    if output_fasta == "":
        raise ValueError("File hasn't got name!")
    
    output_file = output_fasta+".fasta"
    with open(output_file, mode="w") as output_longline:
        output_longline.write(longlines_to_file)
        print(f"FASTA-seq write in { output_file} file")


def read_lines_from_gbk(input_gbk:str, line_name:str) -> list:
    """
    Read lines in GBK and find lines with common pattern and 
    create a list with strings, which start with this pattern. 

    Arguments:
    - input_gbk (str): name of gbk file, where located information about genes
    - line_name (str): pattern, according to which find and add to list lines in GBK

    Example:
    read_lines_from_gbk('example_gbk.gbk', '                     /gene=')[:50] ->
    -> ['phrB', 'dtpD', 'ybgI', 'pxpB', 'pxpC', 'pxpA', 'nei', 'ybgD_1', 'gltA',
        'ybgU', 'sdhC', 'sdhD', 'sdhA', 'sdhB', 'sucA', 'sucB', 'sucC', 'sucD',
        'mngR_1', 'mngA', 'mngB', 'cydA', 'cydB', 'cydX', 'ybgC', 'tolQ', 'tolR',
        'tolA', 'tolB', 'pal_1', 'cpoB', 'nadA', 'pnuC', 'zitB', 'aroG', 'gpmA',
        'galM', 'galK', 'galT', 'galE', 'modF', 'modE', 'acrZ', 'modA', 'modB',
        'btuD_1', 'ybhA', 'pgl', 'oxyR_1', 'ybhH']
    """
    with open(os.path.abspath(input_gbk)) as gbk_file:
        
        list_of_lines = []
        for line in gbk_file.readlines():
            if line.startswith(line_name):
                list_of_lines.append(line.replace(line_name, "").strip("\n").replace('"', ''))

    return list_of_lines

def find_translation_in_gbk(input_gbk:str, gene:str) -> str:
    """
    Find translation-seq for one gene and return it.

    Arguments:
    - input_gbk (str): name of gbk file, where located information about genes
    - gene (str): one gene
    """
    with open(os.path.abspath(input_gbk)) as gbk:
        
        lines = []
        for line in gbk.readlines():
            lines.append(line.replace('\n', '').replace("                     ", " "))
        
        information_about_cds = "".join(lines[18:]).split("     CDS             ")

        cur_gene = '/gene="'+gene+'"'

        for info in information_about_cds:
            if cur_gene in info:

                index_of_translation = info.index('/translation="')

                index_of_last_symbol = 0
                slice_of_info = info[index_of_translation+len('/translation="'):]
                while slice_of_info[index_of_last_symbol] != '"':
                    index_of_last_symbol += 1

                cur_translation = info[index_of_translation+len('/translation="'):index_of_translation+len('/translation="')+index_of_last_symbol]

                return "".join(cur_translation).replace(" ", "")

def find_flanking_genes(input_gbk:str, genes:list) -> list:
    """
    Find translation-seq for each of flanking gene and return a list 
    of translation-seqs.

    Arguments:
    - input_gbk (str): name of gbk file, where located information about genes
    - genes (list): list of a part of flanking genes (before or after)
    """
    list_of_translations = []
    for gene in genes:
        list_of_translations.append(find_translation_in_gbk(input_gbk, gene))

    return list_of_translations


def select_flanking_genes(input_gbk:str, gene:str, n_before:int, n_after:int):
    """
    Select flanking genes for one specific gene and return list of these genes.

     Arguments:
    - input_gbk (str): name of gbk file, where located information about genes
    - genes (list): main genes, for which find flanking genes
    - n_before (int): amount of flanking genes before main gene
    - n_after (int): amount of flanking genes after main gene

    Example:
    select_flanking_genes('example_gbk.gbk', 'cspG_2', 1, 1) -> 
    -> ['MSRKMTGIVKTFDRKSGKGFIIPSDGRKEVQVHISAFTPRDAEVLIPGLRVEFCRVNGLRGPTAANVYLS',
        'msnkmtglvkwfnadkgfgfitpddgskdvfvhftaiqsnefrtlnenqkvefsieqgqrgpaaanvvtl',
        'MNIEELKKQAETEIADFIAQKIAELNKNTGKEVSEIRFTAREKMTGLESYDVKIKIM']
    """
    list_of_genes = read_lines_from_gbk(input_gbk, "                     /gene=")
    
    gene_index = list_of_genes.index(gene)
    genes_before_slice = list_of_genes[gene_index-n_before:gene_index]
    genes_after_slice = list_of_genes[gene_index+1:gene_index+n_after+1]
    
    translation_before = find_flanking_genes(input_gbk, genes_before_slice)
    translation_after = find_flanking_genes(input_gbk, genes_after_slice)
    main_gene = find_translation_in_gbk(input_gbk, gene)

    flanking_genes = translation_before + [main_gene.lower()] + translation_after
    
    return flanking_genes

def select_genes_from_gbk_to_fasta(input_gbk:str, genes:list, n_before=1, n_after=1, output_fasta="") -> str:
    """
    Select flanking genes from GBK for everyone in list of genes (genes). 
    Selected genes write in fasta file.  
    If function works correct, it prints message "Flanking genes write in {output_file} file".
    If output_fasts is empty, raise ValueError("File hasn't got name!").

    Arguments:
    - input_gbk (str): name of gbk file, where located information about genome
    - genes (list): main genes, for which find flanking genes
    - n_before (int): amount of flanking genes before main gene, has default value equal 1
    - n_after (int): amount of flanking genes after main gene, has default value equal 1
    - output_fasta (str): name of file, where saved lines, has default value equal ""

    Structure of result fasta file:
    - every line in upper case is flanking genes of one gene, which write in lower case
    - between groups of flanking genes used empty line
    """    
    if output_fasta == "":
        raise ValueError("File hasn't got name!")

    ans = []
    if type(genes) == list:
        for gene in genes:
            ans.append(select_flanking_genes(input_gbk, gene, n_before, n_after))
    elif type(genes) == str:
        ans.append(select_flanking_genes(input_gbk, genes, n_before, n_after))

    output_file = output_fasta+".fasta"
    with open(output_file, mode="w") as output_fasta:

        for item in ans:
            output_fasta.write("\n".join(item))
            output_fasta.write('\n\n')
        print(f"Flanking genes write in {output_file} file")


 

        

     


