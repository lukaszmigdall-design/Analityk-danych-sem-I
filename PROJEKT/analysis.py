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


