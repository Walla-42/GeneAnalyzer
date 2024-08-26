from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

class DNA:
    def __init__(self, fasta_file=None, seq=None):
        self.sequence = self.parse_sequence(fasta_file, seq)
    
    def parse_sequence(self, File_path=None, seq=None):
        if File_path:
            with open(File_path, "r") as file:
                record = SeqIO.read(file, "fasta")
                return Seq(str(record.seq))
        elif seq:
            return Seq(seq)
        else:
            raise ValueError('Either a file path or a sequence must be provided.')
    
    def GC_count(self):
        gc_count = self.sequence.upper().count('G') + self.sequence.upper().count('C')
        gc_percent = (gc_count / len(self.sequence)) * 100
        return gc_percent, gc_count, len(self.sequence)
    
    def transcribe(self):
        return self.sequence.transcribe()
    
    def reverse_complement(self):
        return self.sequence.reverse_complement()
    
    def base_count(self):
        return Counter(self.sequence.upper())

class RNA:
    def __init__(self, fasta_file=None, seq=None):
        self.sequence = self.parse_sequence(fasta_file, seq)
    
    def parse_sequence(self, File_path=None, seq=None):
        if File_path:
            with open(File_path, "r") as file:
                record = SeqIO.read(file, "fasta")
                return Seq(str(record.seq))
        elif seq:
            return Seq(seq)
        else:
            raise ValueError('Either a file path or a sequence must be provided.')
    
    def GC_count(self):
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        gc_percent = (gc_count / len(self.sequence)) * 100
        return gc_percent
    
    def translate(self):
        return self.sequence.translate()
    
    def base_count(self):
        return Counter(self.sequence.upper())

class Protein:
    def __init__(self, fasta_file=None, seq=None):
        self.sequence = self.parse_sequence(fasta_file, seq)
    
    def parse_sequence(self, File_path=None, seq=None):
        if File_path:
            with open(File_path, "r") as file:
                record = SeqIO.read(file, "fasta")
                return Seq(str(record.seq))
        elif seq:
            return Seq(seq)
        else:
            raise ValueError('Either a file path or a sequence must be provided.')
    
    def protein_to_RNA(self):
        return self.sequence.reverse_translate()
    
    def most_common_amino_acid(self):
        return Counter(self.sequence).most_common(1)[0][0]

if __name__ == "__main__":
    while True:
        type_input = input('What sequence type will you be analyzing (DNA, RNA, Protein)? ').strip().lower()
        
        if type_input not in ['dna', 'rna', 'protein']:
            print('INVALID INPUT TYPE. VALID ENTRIES INCLUDE: DNA, RNA, OR PROTEIN')
            continue 

        file_path = input('Enter file path. Else, type "None": ').strip().lower()

        if file_path == 'none':
            file_path = None

        seq = None
        if file_path is None:
            seq = input('Input your sequence: ')

        if type_input == 'dna':
            sequence = DNA(fasta_file=file_path, seq=seq)
            print(f"GC Content: {sequence.GC_count()}%")
            print(f"Transcribed RNA: {sequence.transcribe()}")
            print(f"Reverse Complement Strand: {sequence.reverse_complement()}")
            print(f"Base Counts of DNA Sequence: {sequence.base_count()}")
        elif type_input == 'rna':
            sequence = RNA(fasta_file=file_path, seq=seq)
            print(f"GC Content: {sequence.GC_count()}%")
            print(f"Translated Protein: {sequence.translate()}")
            print(f"Base Counts of RNA Sequence: {sequence.base_count()}")
        elif type_input == 'protein':
            sequence = Protein(fasta_file=file_path, seq=seq)
            print(f"RNA Sequence: {sequence.protein_to_RNA()}")
            print(f"Most Common Amino Acid: {sequence.most_common_amino_acid()}")

        if input('Would you like to analyze another sequence (yes/no)? ').strip().lower() != 'yes':
            break