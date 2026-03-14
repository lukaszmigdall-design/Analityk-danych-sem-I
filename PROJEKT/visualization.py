import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np


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

def draw_heatmap(frame, df, motifs):

    for widget in frame.winfo_children():
        widget.destroy()

    data = df[motifs].values

    fig, ax = plt.subplots(figsize=(8, 4))

    heatmap = ax.imshow(data, aspect="auto")

    ax.set_xlabel("Segment")
    ax.set_ylabel("Motyw")

    ax.set_yticks(np.arange(len(motifs)))
    ax.set_yticklabels(motifs)

    ax.set_title("Heatmapa wystąpień motywów")

    fig.colorbar(heatmap, ax=ax)

    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill="both", expand=True)

    return canvas

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


def draw_genome_map(frame, df, motifs, islands):

    for widget in frame.winfo_children():
        widget.destroy()

    fig, ax = plt.subplots(figsize=(10, 4))

    segments = df["Segment"]

    # ===== GC content =====
    ax.plot(segments,
            df["GC_content_%"],
            label="GC %",
            linewidth=2)

    # ===== Motywy =====
    for motif in motifs:
        ax.plot(segments,
                df[motif],
                label=motif,
                linestyle="--")

    # ===== CpG islands =====
    for island in islands:

        start_seg = island["Start"] // 100
        end_seg = island["End"] // 100

        ax.axvspan(start_seg,
                   end_seg,
                   alpha=0.3)

    ax.set_xlabel("Segment DNA")
    ax.set_ylabel("Wartość")

    ax.set_title("Mapa genomu")

    ax.legend(loc="upper right")

    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill="both", expand=True)

    return canvas