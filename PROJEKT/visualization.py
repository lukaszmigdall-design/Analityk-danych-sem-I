import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def plot_motif_distribution(df, motif):
    plt.figure()
    plt.bar(df["Segment"], df["Motif_count"])
    plt.xlabel("Segment (po 100 nukleotydów)")
    plt.ylabel("Liczba motywów")
    plt.title(f"Rozmieszczenie motywu '{motif}'")
    plt.tight_layout()
    plt.show()

def draw_plot(root_frame, df, motifs, old_canvas=None):

    if old_canvas:
        old_canvas.get_tk_widget().destroy()

    fig = Figure(figsize=(6, 4))
    ax = fig.add_subplot(111)

    for motif in motifs:
        ax.plot(df["Segment"], df[motif], marker='o', label=motif)

    ax.set_title("Rozmieszczenie motywów")
    ax.set_xlabel("Segment")
    ax.set_ylabel("Liczba wystąpień")
    ax.legend()

    canvas = FigureCanvasTkAgg(fig, master=root_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill="both", expand=True)

    return canvas