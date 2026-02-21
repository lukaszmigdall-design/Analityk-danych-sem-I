# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

#program do analizy motywów DNA
import re
import requests
from Bio import SeqID
from Bio import Entrez
import matplotlib.pyplot as plt
from collections import defaultdict

class DNAmotyw:
    def __init__(self, sequence):
        self.sequece = sequence.upper()
        self.lenght = len(self.sequence)

    @staticmethod
    def from_fasta(file_path):
        record = SeqID.read(file_path, "fasta")
        return DNAmotyw(str(record.seq))

    @staticmethod:

    asdadasdas
    1

    print(

    )

