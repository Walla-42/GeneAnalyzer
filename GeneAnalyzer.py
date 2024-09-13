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
    
    def most_common_amino_acid(self):
        return Counter(self.sequence).most_common(1)[0][0]
    
    def peptide_solubility(self):
        """Function to measure the Grand Average Hydropathy Index  Score (GRAVY index)
        of a peptide of sequence (seq)
        
        input:
        seq: sting
            valid inputs include the single letter codes for the 20 amino acids
            
        return:
        result: int
            average hydropathy index score
            
        NOTE: THIS IS AN AVERAGE AND MAY NOT REFLECT TRUE SOLUBILITY OF THE SEQUENCED PEPTIDE
        """
        sol = {'G': -0.4, 'A': 1.8, 'P': -1.6, 
            'V': 4.2, 'L': 3.8, 'I': 4.5, 
            'M': 1.9, 'F': 2.8, 'Y': -1.3, 
            'W': -0.9, 'S': -0.8, 'T': -0.7,
            'C': 2.5, 'N': -3.5, 'Q': -3.5,
            'K': -3.9, 'H': -3.2, 'R': -4.5,
            'D': -3.5, 'E': -3.5}

        sol_calculate = []
        
        for i in self.sequence:
            if i in sol:
                sol_calculate.append(sol[i])
            else:
                print(f"Unknown amino acid: {i}")
        
        result = sum(sol_calculate)/len(self.sequence)
        
        # print(f'Solubility coefficient values: {sol_calculate}')

        return result


if __name__ == "__main__":
    input_break = '*' * 80
    while True:
        type_input = input('What sequence type will you be analyzing (DNA, RNA, Protein)? ').strip().lower()
        
        if type_input not in ['dna', 'rna', 'protein']:
            print('INVALID INPUT TYPE. VALID ENTRIES INCLUDE: DNA, RNA, OR PROTEIN')
            continue 

        
        while True:
            file_path = input('Enter file path. Else, type "None": ').strip().lower()
            if file_path == 'none':
                file_path = None
                break
            else:
                continue
                

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
            
            while True:
                print(input_break)
                print('RNA Menu:\n1. GC Content\n2. Translate RNA Sequence\n3. Base Count Analysis')
                rna_menu_input = input('What type of analysis would you like to perform (select 1-4): ')
                print(input_break)
                if rna_menu_input.strip() == '1':
                     print(f"GC Content: {rna.GC_count()}%")
                elif rna_menu_input.strip() == '2':
                    print(f"Translated Protein: {rna.translate()}")
                elif rna_menu_input.strip() == '3':
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
            while True:
                print(input_break)
                print('Protein Analysis Menu\n1. Estimate Solubility of Peptide\n2. Most Common Amino Acid')
                protein_menu_input = input('What type of analysis do you want to perform (select 1-2): ')
                print(input_break)
                if protein_menu_input.strip() == '1': 
                    result = protein.peptide_solubility()
                    if result >= 0:
                        print(f'Peptide is soluble in water with a grand average of hydropathy index (GRAVY) score of {round(result,3)}')
                    else:
                        print(f'Peptide is not soluble in water with a grand average of hydropathy index (GRAVY) score of {round(result,3)}')
                elif protein_menu_input.strip() == '2':
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
