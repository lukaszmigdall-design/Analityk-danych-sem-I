import numpy as np
import pandas as pd


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