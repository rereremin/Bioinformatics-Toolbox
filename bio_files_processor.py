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
