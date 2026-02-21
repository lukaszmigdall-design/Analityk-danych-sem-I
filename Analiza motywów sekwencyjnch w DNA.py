# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

#program do analizy motywów DNA
from Bio import SeqID, SeqIO
from Bio import Entrez


class DNAmotyw:
    def __init__(self, sequence):
        self.sequece = sequence.upper()
        self.lenght = len(self.sequence)

    @staticmethod
    def from_fasta(file_path):
        record = SeqID.read(file_path, "fasta")
        return DNAmotyw(str(record.seq))

    @staticmethod
    def from_ncbi(accesion, email):
        Entrez.email = email
        handle = Entrez.efetch(db="nucleotide", id=accesion, rettype='fasta', retmode='text')
        record = SeqIO.read(handle, "fasta")
        record.close()
        return DNAmotyw(str(record.seq))

def main():
    print(Analiza motywów sekwencjnych w DNA)
    choice = input("Opcja 1 - wczytaj z pliku FASTA\nOpcja 2 - wczytaj z NCBI\nWybór(1 lub 2): ")
    if choice == "1":
        file_path =



