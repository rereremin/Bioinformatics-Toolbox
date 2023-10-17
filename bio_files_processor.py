import os

def convert_multiline_fasta_to_oneline(input_fasta:str, output_fasta="") -> str:
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
        print(f"FASTA-seq write in {output_file} file")

def read_lines_from_gbk(input_gbk:str, line_name:str) -> list:
    """
    
    """
    with open(os.path.abspath(input_gbk)) as gbk_file:
        
        list_of_lines = []
        for line in gbk_file.readlines():
            if line.startswith(line_name):
                list_of_lines.append(line.replace(line_name, "").strip("\n").replace('"', ''))

    return list_of_lines

def find_translation_in_gbk(input_gbk:str, gene:str):
    """
    
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

def find_flanking_genes(input_gbk:str, genes:str) -> list:
    """
    
    """
    list_of_translations = []
    for gene in genes:
        list_of_translations.append(find_translation_in_gbk(input_gbk, gene))

    return list_of_translations


def select_flanking_genes(input_gbk:str, gene:str, n_before:int, n_after:int):
    """
    
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
