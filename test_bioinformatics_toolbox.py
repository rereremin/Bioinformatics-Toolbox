import pytest
from bio_files_processor import convert_multiline_fasta_to_oneline
from bioinformatics_toolbox import DNASequence, RNASequence, AminoAcidSequence


class TestDNASequence:
    """
    Test methods of the DNASequence class.
    """
    def test_transcribe(self):
        """
        Test the transcribe method of the DNASequence class with mixed case input.
        """
        dna_seq = DNASequence("ATTTTCGCtttcgcGCT", "DNA")
        assert dna_seq.transcribe() == "AUUUUCGCuuucgcGCU"

    def test_gc_content(self):
        """
        Test the GC content method of the DNASequence class with mixed case input.
        """
        dna_seq = DNASequence("ATTTTCGCtttcgcGCTGGGGCCGC", "DNA")
        gc_content = dna_seq.gc_content()
        assert gc_content == 52.0


class TestRNASequence:
    """
    Test methods of the RNASequence class.
    """
    def test_translate_stop_codon(self):
        """
        Test raising ValueError for an invalid RNA 
        sequence with incorrect stop codon.
        """
        rna_seq = RNASequence("AUGUAAUAGCCC", "RNA")
        with pytest.raises(ValueError):
            rna_seq.translation()
    
    def test_translate_triplets_num(self):
        """
        Test raising ValueError for an invalid RNA sequence with incorrect length.
        """
        rna_seq = RNASequence("AUGGGGCGCG", "RNA")
        with pytest.raises(ValueError):
            rna_seq.translation()

class TestAminoAcidSequence:
    """
    Test methods of the AminoAcidSequence class.
    """
    def test_compute_hydrophobicity(self):
        """
        Test the compute hydrophobicity method of the AminoAcidSequence class.
        """
        aminoacid_seq = AminoAcidSequence("ACHTTTTVLI", "protein")
        hydrophobicity = aminoacid_seq.compute_hydrophobicity()
        assert hydrophobicity == 40.0

    def test_compute_(self):
        """
        Check the type method of the AminoAcidSequence class.
        """
        aminoacid_seq = AminoAcidSequence("ACHTTTTVLI", "protein")
        protein_type = aminoacid_seq.is_type()
        assert protein_type == True
    


class TestFastaProcessing:
    """
    Test functions related to processing FASTA files.
    """

    @pytest.fixture
    def multiline_fasta_content(self) -> str:
        """
        Fixture providing multiline FASTA content.
        """

        return ">Sequence1 Species1\nATGC\nCAATCG\nGAT\n>Sequence2 Species2\nTTAA\nCCGG\n>Sequence3 Species3\nGATTACA\n"

    @pytest.fixture
    def multiline_fasta_file(self, tmp_path: str, multiline_fasta_content: str) -> str:
        """
        Fixture creating a temporary multiline FASTA file.
        """
        filepath = tmp_path / "multiline.fasta"
        with open(filepath, "w") as f:
            f.write(multiline_fasta_content)
        return filepath

    def test_convert_multiline_fasta_to_oneline_output(self, multiline_fasta_file: str, tmp_path: str) -> None:
        """
        Test conversion of multiline FASTA to oneline format and check output file existence.
        """
        output_file = tmp_path / "one_line.fasta"
        convert_multiline_fasta_to_oneline(multiline_fasta_file, str(output_file))
        assert output_file.exists()

    def test_convert_multiline_fasta_to_oneline_content(self, multiline_fasta_file: str, tmp_path: str) -> None:
        """
        Test conversion of multiline FASTA to oneline format and check output content.
        """
        output_file = tmp_path / "one_line.fasta"
        convert_multiline_fasta_to_oneline(multiline_fasta_file, str(output_file))
        with open(output_file, "r") as f:
            output_content = f.read()
            assert output_content == ">Sequence1 Species1\nATGCCAATCGGAT\n>Sequence2 Species2\nTTAACCGG\n>Sequence3 Species3\nGATTACA\n"

