import numpy as np
import pandas as pd




#poprawność danych tj czy mamy DNA
def validate_dna(sequence):
    allowed = set("ATGC")
    sequence = sequence.upper()

    for base in sequence:
        if base not in allowed:
            return False
    return True


#motyw
def find_motif(sequence, motif):
    positions = []
    m_len = len(motif)

    for i in range(len(sequence) - m_len + 1):
        if sequence[i:i+m_len] == motif:
            positions.append(i)

    return positions

#GC-content
def gc_content(sequence):
    g = sequence.count('G')
    c = sequence.count('C')
    total = len(sequence)

    if total == 0:
        return 0
    return ((g + c) / total * 100, 2)

def segment_sequence_multiple(sequence, motifs, segment_size=100):

    segments = [sequence[i:i+segment_size]
                for i in range(0, len(sequence), segment_size)]

    data = {
        "Segment": np.arange(len(segments)),
        "Start": np.arange(0, len(sequence), segment_size)
    }


#motywy
    for motif in motifs:
        counts = [seg.count(motif) for seg in segments]
        data[motif] = counts


# GC-content
    gc_values = []
    for seg in segments:
        g = seg.count("G")
        c = seg.count("C")
        total = len(seg)

        if total == 0:
            gc = 0
        else:
            gc = (g + c) / total * 100

        gc_values.append(round(gc, 2))

    data["GC_content_%"] = gc_values

    df = pd.DataFrame(data)
    return df

#wyspy CpG
def find_cpg_islands(sequence, window_size=200):

    islands = []

    for i in range(len(sequence) - window_size):

        window = sequence[i:i+window_size]

        c = window.count("C")
        g = window.count("G")
        cg = window.count("CG")

        if len(window) == 0:
            continue

        gc_content = (c + g) / len(window)

        expected = (c * g) / len(window) if len(window) > 0 else 0

        if expected == 0:
            ratio = 0
        else:
            ratio = cg / expected

        if gc_content >= 0.5 and ratio >= 0.6:

            islands.append({
                "Start": i,
                "End": i + window_size,
                "GC_content": round(gc_content * 100, 2),
                "CpG_ratio": round(ratio, 2)
            })

    return islands


#Transkrypcja i translacja
def dna_to_rna(sequence):
    sequence = sequence.upper()
    rna = sequence.replace("T", "U")
    return rna

codon_table = {

"AUG":"M",

"UUU":"F","UUC":"F",
"UUA":"L","UUG":"L",
"CUU":"L","CUC":"L","CUA":"L","CUG":"L",

"AUU":"I","AUC":"I","AUA":"I",

"GUU":"V","GUC":"V","GUA":"V","GUG":"V",

"UCU":"S","UCC":"S","UCA":"S","UCG":"S",

"CCU":"P","CCC":"P","CCA":"P","CCG":"P",

"ACU":"T","ACC":"T","ACA":"T","ACG":"T",

"GCU":"A","GCC":"A","GCA":"A","GCG":"A",

"UAU":"Y","UAC":"Y",

"CAU":"H","CAC":"H",

"CAA":"Q","CAG":"Q",

"AAU":"N","AAC":"N",

"AAA":"K","AAG":"K",

"GAU":"D","GAC":"D",

"GAA":"E","GAG":"E",

"UGU":"C","UGC":"C",

"UGG":"W",

"CGU":"R","CGC":"R","CGA":"R","CGG":"R",

"AGU":"S","AGC":"S",

"AGA":"R","AGG":"R",

"GGU":"G","GGC":"G","GGA":"G","GGG":"G",

"UAA":"*","UAG":"*","UGA":"*"
}

def rna_to_protein(rna):
    protein = ""
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        amino_acid = codon_table.get(codon, "X")
        protein += amino_acid
    return protein


def find_orfs(rna):

    start = "AUG"
    stops = ["UAA","UAG","UGA"]

    orfs = []

    for i in range(0,len(rna)-3):

        codon = rna[i:i+3]

        if codon == start:

            for j in range(i+3,len(rna)-3,3):

                stop = rna[j:j+3]

                if stop in stops:

                    orfs.append((i,j+3))
                    break

    return orfs