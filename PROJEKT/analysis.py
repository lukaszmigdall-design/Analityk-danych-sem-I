import numpy as np
import pandas as pd

def find_motif(sequence, motif):
    positions = []
    m_len = len(motif)

    for i in range(len(sequence) - m_len + 1):
        if sequence[i:i+m_len] == motif:
            positions.append(i)

    return positions


def segment_sequence(sequence, motif, segment_size=100):
    segments = [sequence[i:i+segment_size]
                for i in range(0, len(sequence), segment_size)]

    counts = []
    for seg in segments:
        count = 0
        for i in range(len(seg) - len(motif) + 1):
            if seg[i:i+len(motif)] == motif:
                count += 1
        counts.append(count)

    df = pd.DataFrame({
        "Segment": np.arange(len(segments)),
        "Start": np.arange(0, len(sequence), segment_size),
        "Motif_count": counts
    })

    return df