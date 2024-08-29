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
        return round(gc_percent, 2)
    
    def transcribe(self):
        if set('ATCG').issubset(set(self.sequence)) is False:
            print('Warning: Some bases not recognized and have not been transcribed')
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
        gc_count = self.sequence.upper().count('G') + self.sequence.upper().count('C')
        gc_percent = (gc_count / len(self.sequence)) * 100
        return round(gc_percent, 2)
    
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
    input_break = '*' * 80
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
            dna = DNA(fasta_file=file_path, seq=seq)
            while True:
                print(input_break)
                print('DNA Menu:\n1. GC Content\n2. Transcibe DNA\n3. Reverse Compliment Strand\n4. Base Count Analysis')
                dna_menu_input = input('What type of analysis would you like to perform (select 1-4): ')
                print(input_break)
                if dna_menu_input.strip() == '1':
                    print(f"GC Content: {dna.GC_count()}%")
                elif dna_menu_input.strip() == '2':
                    print(f"Transcribed RNA: {dna.transcribe()}")
                elif dna_menu_input.strip() == '3':
                    print(f"Reverse Complement Strand: {dna.reverse_complement()}")
                elif dna_menu_input.strip() == '4':
                    print(f"Base Counts of DNA Sequence: {dna.base_count()}")
                else:
                    print('INVALID INPUT TYPE. VALID ENTRIES INCLUDE: 1, 2, 3, OR 4')
                    print(input_break)
                    continue
                print(input_break)
                if input('Would you like to perform another analysis on this sequence? (yes/no) ').strip().lower() == 'no':
                    break

        elif type_input == 'rna':
            rna = RNA(fasta_file=file_path, seq=seq)
            print(input_break)
            print('RNA Menu:\n1. GC Content\n2. Translate RNA Sequence\n3. Base Count Analysis')
            rna_menu_input = input('What type of analysis would you like to perform (select 1-4): ')
            print(input_break)
            while True:
                if rna_menu_input.strip() == 1:
                     print(f"GC Content: {rna.GC_count()}%")
                elif rna_menu_input.strip() == 2:
                    print(f"Translated Protein: {rna.translate()}")
                elif rna_menu_input.strip() == 3:
                    print(f"Base Counts of RNA Sequence: {rna.base_count()}")
                else:
                    print('INVALID INPUT TYPE. VALID ENTRIES INCLUDE: 1, 2, or 3')
                    print(input_break)
                    continue

                print(input_break)
                if input('Would you like to perform another analysis on this sequence? (yes/no) ').strip().lower() == 'no':
                    break
            
        elif type_input == 'protein':
            protein = Protein(fasta_file=file_path, seq=seq)
            print(input_break)
            print('Protein Analysis Menu\n1. Reverse Translation\n2. Most Common Amino Acid')
            protein_menu_input = input('What type of analysis do you want to perform (select 1-2): ')

            print(input_break)
            if protein_menu_input.strip() == 1: 
                print(f"RNA Sequence: {protein.protein_to_RNA()}")
            elif protein_menu_input.strip() == 2:
                print(f"Most Common Amino Acid: {protein.most_common_amino_acid()}")
            else:
                    print('INVALID INPUT TYPE. VALID ENTRIES INCLUDE: 1, or 2')
                    print(input_break)
                    continue

            print(input_break)
            if input('Would you like to perform another analysis on this sequence? (yes/no) ').strip().lower() == 'no':
                break
            
        if input('Would you like to analyze another sequence (yes/no)? ').strip().lower() != 'yes':
            break
