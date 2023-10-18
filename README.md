# Bioinformatics-Toolbox

**Bioinformatics-Toolbox** - is a tool which allows the performing diffrent procedures with DNA/RNA, protein and filter fastq-libraries. 

### Usage

The tool works by calling the functions `run_dna_rna_tools`,`run_protein_tools` which takes arbitrary number of arguments with nucleic acids and protein sequencies (*str*) and the name of the procedure to be performed (always the last argument, *str*, see the usage examples below). The output is the result of the procedure as *string or tuple* if one sequence is submitted or *list* if several. Tool include`filter_fastq`, which takes fastq-dict and three characteristics to filtration (*dict*, *tuple*, *tuple*, *int*) .The output is the result of the dictionary filtration as *dict*. Also tool include `convert_multiline_fasta_to_oneline` and `select_genes_from_gbk_to_fasta`, these functions able to work with fasta and gbk (first function takes multiline fasta file, name of output file and returns oneline fasta file; second function takes GBK file, genes, and amount of genes before and after, output file name and returns fasta file with flanking genes groups).

**NOTE:**  For the procedure `check_mutations` in `run_protein_tools` a fixed number of string arguments are used: one RNA sequence, one protein sequence and the name of procedure itself.

### Modules
#### dna_rna_tools
- `reverse_complement` - perform reverse complementary complementary sequence of RNA or DNA
- `complement_seq` - perform complementary complementary sequence of RNA or DNA.
- `reverse_seq` - perform reverse DNA or RNA sequence
- `transcribe_dna` - perform the transcribtion of DNA into mRNA

#### protein_tools
- `compute_molecular_weight` — computes molecular weight of protein sequence in g/mol
- `compute_length` — computes the number of amino acids in protein sequence
- `compute_hydrophobicity` — computes the percentage of gydrophobic aminoacids in protein sequence
- `check_mutations` — checks missense mutations in the protein sequence after translation
- `protein_to_dna`- returns possible variants of DNAs for a given protein sequence
- `count_amino_acids` - calculates the number of each aminoacid in protein sequence

#### fastq_filtration
- `filter_by_quality` -  filters fastq-sequence by quality and decoding quality sequence into q-score list and calculate quality
- `filtre_by_length` - filters fastq-sequence by length
- `filter_by_gc_level` - filters sequence by gc-level
- `save_filtered_fastq` - Saves filtered fastq-sequencies and their characteristics in directory, which called fastq_filtrator_resuls
- `read_fastq_file` - reads lines in fastq-file

  **NOTE:** This five functions use together to filter fastq-library and write it in file.

#### bio_files_processor
- `convert_multiline_fasta_to_oneline` - converts multiline fasta to oneline and writes in fasta file
- `read_lines_from_gbk` - reads lines in GBK and find lines with common pattern and creates a list with strings, which start with this pattern
- `find_translation_in_gbk` - find translation-seq for one gene and return it
- `find_flanking_genes` - find translation-seq for each of flanking gene and return a list of translation-seqs
- `select_flanking_genes` - select flanking genes for one specific gene and return list of these genes
  
  **NOTE:** `read_lines_from_gbk`, `find_translation_in_gbk`, `find_flanking_genes`, `select_flanking_genes` supporting functions for `select_genes_from_gbk_to_fasta`
  

### Examples
#### dna_rna_tools
```python
run_dna_rna_tools("ATCATAGA","CAGatAGC","transcribe_dna")
# ['AUCAUAGA', 'CAGauAGC']
run_dna_rna_tools("UAVGCGCGC","transcribe_dna")
# ValueError: UAVGCGCGC isn't RNA or DNA.
run_dna_rna_tools("UAGCGCGC","transcribe_dna")
# ValueError: UAVGCGCGC isn't RNA or DNA.

run_dna_rna_tools("ATCATAGA","CAGatAGC","UUUAUAUA", "reverse")
# # ['AGATACTA', 'CGAtaGAC', 'AUAUAUUU']
run_dna_rna_tools("ACACARRR", "reverse")
# ValueError: ACACARRR isn't RNA or DNA.

run_dna_rna_tools("ATCATAGA","CAGatAGC","UUUAUAUA", "complement")
# # ['TAGTATCT', 'GTCtaTCG', 'AAAUAUAU']
run_dna_rna_tools("ACACARRR", "complement")
# # ValueError: ACACARRR isn't RNA or DNA.

run_dna_rna_tools("ATCATAGA","CAGatAGC","UUuAUAUA", 'reverse_complement')
# # ['TCTATGAT', 'GCTatCTG', 'UAUAUaAA']
run_dna_rna_tools("ACACARRR", "reverse_complement")
# #  ValueError: ACACARRR isn't RNA or DNA.
```
#### protein_tools
```python
run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_length')
#[('MAEGEITNLP', 10), ('tGQYLAMDTSgLLYGSQT', 18)]

run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_molecular_weight')
#[('MAEGEITNLP', 1055.496), ('tGQYLAMDTSgLLYGSQT', 1886.872)]

run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_hydrophobicity')
#[('MAEGEITNLP', 50.0), ('tGQYLAMDTSgLLYGSQT', 27.778)]

run_protein_tools('AUGGAUCAUcAAUAA', 'MDKL*', 'check_mutations')
#'Mutations: K3, L4.'

run_protein_tools('MAEGLP', 'LYGSQT','protein_to_dna')
#['ATG GCT/GCC/GCA/GCG GAA/GAG GGT/GGC/GGA/GGG TTA/TTG/CTT/CTC/CTA/CTG CCT/CCC/CCA/CCG',
#'TTA/TTG/CTT/CTC/CTA/CTG TAT/TAC GGT/GGC/GGA/GGG TCT/TCC/TCA/TCG/AGT/AGC CAA/CAG ACT/ACC/ACA/ACG']

run_protein_tools('MAEGLP', 'LYGSQT','count_amino_acids')
#[{'M': 1, 'A': 1, 'E': 1, 'G': 1, 'L': 1, 'P': 1},
#{'L': 1, 'Y': 1, 'G': 1, 'S': 1, 'Q': 1, 'T': 1}]
```
#### fastq_filtration
```python
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
d = filter_fastq(EXAMPLE_FASTQ, (20, 80), (0, 80), 30)
print(d)

#{'@SRX079804:1:SRR292678:1:1101:170868:170868': ('CTGCCGAGACTGTTCTCAGACATGGAAAGCTCGATTCGCATACACTCGCTGAGTAAGAGAGTCACACCAAATCACAGATT', #'E;FFFEGFGIGGFBG;C6D<@C7CDGFEFGFHDFEHHHBBHHFDFEFBAEEEEDE@A2=DA:??C3<BCA7@DCDEG*EB'),
#'@SRX079804:1:SRR292678:1:1101:171075:171075': ('CATTATAGTAATACGGAAGATGACTTGCTGTTATCATTACAGCTCCATCGCATGAATAATTCTCTAATATAGTTGTCAT', #'HGHHHHGFHHHHFHHEHHHHFGEHFGFGGGHHEEGHHEEHBHHFGDDECEGGGEFGF<FGGIIGEBGDFFFGFFGGFGF'),
#'@SRX079804:1:SRR292678:1:1101:175500:175500': ('GACGCCGTGGCTGCACTATTTGAGGCACCTGTCCTCGAAGGGAAGTTCATCTCGACGCGTGTCACTATGACATGAATG', #'GGGGGFFCFEEEFFDGFBGGGA5DG@5DDCBDDE=GFADDFF5BE49<<<BDD?CE<A<8:59;@C.C9CECBAC=DE'),
#'@SRX079804:1:SRR292678:1:1101:190136:190136': ('GAACCTTCTTTAATTTATCTAGAGCCCAAATTTTAGTCAATCTATCAACTAAAATACCTACTGCTACTACAAGTATT', #'DACD@BEECEDE.BEDDDDD,>:@>EEBEEHEFEHHFFHH?FGBGFBBD77B;;C?FFFFGGFED.BBABBG@DBBE'),
#'@SRX079804:1:SRR292678:1:1101:190845:190845': ('CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC', #'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE'),
#'@SRX079804:1:SRR292678:1:1101:204480:204480': ('AGTGAGACACCCCTGAACATTCCTAGTAAGACATCTTTGAATATTACTAGTTAGCCACACTTTAAAATGACCCG','<98;<@@@:@CD@BCCDD=DBBCEBBAAA@9???@BCDBCGF=GEGDFGDBEEEEEFFFF=EDEE=DCD@@BBC')}
```
#### convert_multiline_fasta_to_oneline
```python
  convert_multiline_fasta_to_oneline('example_multiline_fasta.fasta', "fasta1")
```
   ([example_fasta.fasta](https://github.com/rereremin/Bioinformatics-Toolbox/blob/hw6/examples/example_multiline_fasta.fasta)) ->  ([fasta1.fasta](https://github.com/rereremin/Bioinformatics-Toolbox/blob/hw6/examples/fasta1.fasta))

#### select_genes_from_gbk_to_fasta
```python
   select_genes_from_gbk_to_fasta('example_gbk.gbk', 'uvrB', 2, 2, 'gbk1')
```
   ([example_gbk.gbk](https://github.com/rereremin/Bioinformatics-Toolbox/blob/hw6/examples/example_gbk.gbk)) -> ([gbk1.fasta](https://github.com/rereremin/Bioinformatics-Toolbox/blob/hw6/examples/gbk1.fasta))

**NOTE:** 
Structure of result fasta file:
- every line in upper case is flanking genes of one gene, which write in lower case
- between groups of flanking genes used empty line
   
### Additional information
- The `dna_rna_tools` works **only** with DNA and RNA sequences. If any of the entered sequences contains inappropriate characters or cannot exist, the program prints message "It isn't RNA or DNA".
- The `protein_tools` works **only** with protein and RNA sequences. If any of the entered sequences contain inappropriate characters or cannot exist, the program will display an error. Sequences can contain characters of any case.
- For the procedure `check_mutations` there are extra requirements for RNA and protein sequences: mRNA sequences must contain **start-codon** and **one of the stop-codons**, protein sequnces must start with **"M"** and ends with **"*"** (stop-codon). 
```python
run_protein_tools("AUGGUAGGGAAAUUUUGA", "MGGKF", 'check_mutations') #ValueError: Stop (*) is absent
run_protein_tools("AUGGUAGGGAAAUUUUGA", "GGKF*", 'check_mutations') #ValueError: Start (M) is absent
```

