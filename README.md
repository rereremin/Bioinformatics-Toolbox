# BreadcrumbsBioinformatics-Toolbox
This repo was created over the course of one academic year while studying at the Bioinformatics Institute. It provides basic functionality for working with **biological sequences**, processing information from **fasta** and **gbk** files, and presents the implementation of a **random forest using parallelisation**. 

Also here you can find tests for checking modules and examples of using functions from modules.

Let's take a detailed look at each of the files presented in this repository.

## bioinformatics_toolsbox.py
This module is written using OOP and contains classes that store methods and attributes of three known biological sequences (DNA, RNA, proteins). 

The structure is represented by class `BiologicalSequence`, from which class `NucleicAcidSequence` is inherited. Class `NucleicAcidSequence` is the parent of classes `DNASequence` and `RNASequence`. And class `AminoAcidSequence` is a child of class `BiologicalSequence`. 

In addition to the classes described above, other functions `filter_fastq` and `telegram_logger` (decorator) are implemented in this module. They are described below: 

- `filter_fastq` filters Fastq format files by three attributes: length, GC-composition, quality. Output filtered sequences to a new file.
- `telegram_logger` this function is implemented on Telegram bot API and it is needed to send messages and log files of running scripts to Telegram chat for notification.
