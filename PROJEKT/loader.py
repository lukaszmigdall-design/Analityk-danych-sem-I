from Bio import SeqIO, Entrez
from tkinter import messagebox

Entrez.email = "twoj_email@example.com"

def load_sequence_from_file(file_path):
    if file_path.endswith((".fasta", ".fa")):
        records = list(SeqIO.parse(file_path, "fasta"))
        sequence = "".join(str(rec.seq) for rec in records)
        return sequence.upper()
    else:
        with open(file_path, "r") as f:
            return f.read().replace("\n", "").upper()


def load_sequence_from_ncbi(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession,
                               rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq).upper()
    except Exception as e:
        messagebox.showerror("Błąd NCBI", str(e))
        return None