# Bioinformatics-Toolbox
This repo was created over the course of one academic year while studying at the Bioinformatics Institute. It provides basic functionality for working with **biological sequences**, processing information from **fasta** and **gbk** files, and presents the implementation of a **random forest using parallelisation**. 

Also here you can find tests for checking modules and examples of using functions from modules.

Let's take a detailed look at each of the files presented in this repository.

## bioinformatics_toolsbox.py
This module is written using OOP and contains classes that store methods and attributes of three known biological sequences (DNA, RNA, proteins). 

The structure is represented by class `BiologicalSequence`, from which class `NucleicAcidSequence` is inherited. Class `NucleicAcidSequence` is the parent of classes `DNASequence` and `RNASequence`. And class `AminoAcidSequence` is a child of class `BiologicalSequence`. 

In addition to the classes described above, other functions `filter_fastq` and `telegram_logger` (decorator) are implemented in this module. They are described below: 

- **func**  `filter_fastq` filters Fastq format files by three attributes: length, GC-composition, quality. Output filtered sequences to a new file.
- **func**  `telegram_logger` this function is implemented on Telegram bot API and it is needed to send messages and log files of running scripts to Telegram chat for notification.

## bio_files_processor.py
This module implements functionality for working with fasta and gbk files:
1) functions and classes, to work with fasta:
   - **func** `convert_multiline_fasta_to_oneline` converts multiline fasta to oneline and writes in fasta file;
   - **class** `FastaRecord` returns a single record in a FASTA file;
   - **class** `OpenFasta` an iterator over fasta files, displaying the id description of the sequence and the sequence itself.
2) functions, to work with gbk:
   - func `read_lines_from_gbk` read lines in GBK and find lines with common pattern and create a list with strings, which start with this pattern;
   - **func** `find_translation_in_gbk` find translation-seq for one gene and return it;
   - **func** `find_flanking_genes` find translation-seq for each of flanking gene and return a list 
    of translation-seqs;
   - **func** `select_flanking_genes` select flanking genes for one specific gene and return list of these genes;
   - **func** `select_genes_from_gbk_to_fasta` select flanking genes from GBK for everyone in list of genes (genes) and write in fasta file.
  
## custom_random_forest.py
This module implements a random forest using parallelisation. 

## test_bioinformatics_toolbox.py
Here are 8 tests that test functions from **bioinformatics_toolsbox.py** and **bio_files_processor.py**. There are 4 types classes with tests: 
1) **class** `TestDNASequence` (test_transcribe, test_gc_content)
2) **class** `TestRNASequence` (test_translate_stop_codon, test_translate_triplets_num)
3) **class** `TestAminoAcidSequence` (test_compute_hydrophobicity, test_type)
4) **class** `TestFastaProcessing` (test_convert_multiline_fasta_to_oneline_output, test_convert_multiline_fasta_to_oneline_content)

## Showcases.ipynb
Examples of using functions from the above modules and running tests.

