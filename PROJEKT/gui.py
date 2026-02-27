import tkinter as tk
from tkinter import filedialog, messagebox

from loader import load_sequence_from_file, load_sequence_from_ncbi
from analysis import find_motif, segment_sequence
from visualization import plot_motif_distribution


class DNAApp:

    def __init__(self, root):
        self.root = root
        self.root.title("Projekt 1 - Analiza motywów DNA")
        self.root.geometry("1500x950")

        self.results_df = None

        tk.Label(root, text="Plik FASTA/TXT:").pack()
        self.file_entry = tk.Entry(root, width=50)
        self.file_entry.pack()

        tk.Button(root, text="Wybierz plik",
                  command=self.browse_file).pack(pady=5)

        tk.Label(root, text="LUB podaj Accession ID (NCBI):").pack()
        self.accession_entry = tk.Entry(root, width=30)
        self.accession_entry.pack()

        tk.Label(root, text="Motyw (np. ATG, TATA):").pack(pady=5)
        self.motif_entry = tk.Entry(root, width=20)
        self.motif_entry.pack()

        tk.Button(root, text="Analizuj",
                  command=self.analyze).pack(pady=10)

        self.result_label = tk.Label(root, text="")
        self.result_label.pack()

        tk.Button(root, text="Eksportuj do CSV",
                  command=self.export_csv).pack(pady=10)

    def browse_file(self):
        filename = filedialog.askopenfilename(
            filetypes=[("FASTA files", "*.fasta *.fa"),
                       ("Text files", "*.txt")]
        )
        self.file_entry.delete(0, tk.END)
        self.file_entry.insert(0, filename)

    def analyze(self):
        motif = self.motif_entry.get().upper()
        file_path = self.file_entry.get()
        accession = self.accession_entry.get()

        if not motif:
            messagebox.showerror("Błąd", "Podaj motyw!")
            return

        if file_path:
            sequence = load_sequence_from_file(file_path)
        elif accession:
            sequence = load_sequence_from_ncbi(accession)
            if sequence is None:
                return
        else:
            messagebox.showerror("Błąd", "Wybierz plik lub podaj Accession ID!")
            return

        positions = find_motif(sequence, motif)
        count = len(positions)

        self.result_label.config(
            text=f"Długość sekwencji: {len(sequence)}\n"
                 f"Liczba wystąpień: {count}"
        )

        self.results_df = segment_sequence(sequence, motif)

        plot_motif_distribution(self.results_df, motif)

    def export_csv(self):
        if self.results_df is None:
            messagebox.showerror("Błąd", "Najpierw wykonaj analizę!")
            return

        save_path = filedialog.asksaveasfilename(
            defaultextension=".csv"
        )

        if save_path:
            self.results_df.to_csv(save_path, index=False)
            messagebox.showinfo("Sukces", "Dane zapisane do CSV")