import matplotlib.pyplot as plt

def plot_motif_distribution(df, motif):
    plt.figure()
    plt.bar(df["Segment"], df["Motif_count"])
    plt.xlabel("Segment (po 100 nukleotydów)")
    plt.ylabel("Liczba motywów")
    plt.title(f"Rozmieszczenie motywu '{motif}'")
    plt.tight_layout()
    plt.show()