import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from analysis import validate_dna
from loader import load_sequence_from_file, load_sequence_from_ncbi
from analysis import find_motif, segment_sequence_multiple, gc_content
from visualization import draw_plot
from visualization import draw_heatmap
from analysis import find_cpg_islands
from visualization import draw_genome_map

class DNAApp:

    def __init__(self, root):

        self.root = root
        self.setup_style()

        self.root.title("Analiza motywów DNA")
        self.root.geometry("1200x900")

        self.sequence = ""
        self.results_df = None
        self.canvas = None

        self.create_widgets()

    def setup_style(self):

        style = ttk.Style()

        # Motyw
        style.theme_use("clam")

        # Globalna czcionka
        style.configure(".",
                        font=("Segoe UI", 10))

        # ===== Zielony przycisk Analizuj =====
        style.configure("Success.TButton",
                        background="#2ecc71",
                        foreground="white",
                        padding=6)

        style.map("Success.TButton",
                  background=[("active", "#27ae60")])

        # ===== Niebieski przycisk Wybierz =====
        style.configure("Primary.TButton",
                        background="#3498db",
                        foreground="white",
                        padding=6)

        style.map("Primary.TButton",
                  background=[("active", "#2980b9")])

        # ===== Pomarańczowy Eksport =====
        style.configure("Warning.TButton",
                        background="#f39c12",
                        foreground="white",
                        padding=6)

        style.map("Warning.TButton",
                  background=[("active", "#d35400")])
        # Przyciski
        style.configure("TButton",
                        padding=6)

        # LabelFrame
        style.configure("TLabelframe",
                        background="#f4f6f9")

        # Treeview
        style.configure("Treeview",
                        rowheight=25)

        # Kolor tła okna
        self.root.configure(bg="#f4f6f9")


    def toggle_motif_list(self):

        if self.motif_visible:
            self.motif_frame.grid_remove()
            self.toggle_button.config(text="Wybierz motywy ▼")
            self.motif_visible = False
        else:
            self.motif_frame.grid(row=3, column=1, sticky="w")
            self.toggle_button.config(text="Ukryj motywy ▲")
            self.motif_visible = True

    def update_motif_list(self, motifs):

        self.motif_listbox.delete(0, tk.END)
        for motif in motifs:
            self.motif_listbox.insert(tk.END, motif)

    def filter_motifs(self, *args):

        search_text = self.search_var.get().upper()

        if search_text == "":
            filtered = self.all_motifs
        else:
            filtered = [m for m in self.all_motifs
                        if search_text in m]

        self.update_motif_list(filtered)

    # =============================
    def create_widgets(self):
            # =========================
            # GŁÓWNY PODZIAŁ OKNA
            # =========================

            main_frame = ttk.Frame(self.root)
            main_frame.pack(fill="both", expand=True)

            # ===== LEWA STRONA (motywy) =====
            left_frame = ttk.Frame(main_frame, width=220)
            left_frame.pack(side="left", fill="y")
            left_frame.pack_propagate(False)

            motif_frame = ttk.LabelFrame(left_frame, text="Motywy DNA")
            motif_frame.pack(fill="both", expand=True, padx=10, pady=10)

            # ===== Pole wyszukiwania =====
            ttk.Label(motif_frame, text="Szukaj:").pack(pady=(5, 0))

            self.search_var = tk.StringVar()
            self.search_var.trace("w", self.filter_motifs)

            self.search_entry = ttk.Entry(motif_frame,
                                          textvariable=self.search_var)
            self.search_entry.pack(fill="x", padx=5, pady=5)

            # lista motywów
            scrollbar = ttk.Scrollbar(motif_frame)
            scrollbar.pack(side="right", fill="y")

            self.motif_listbox = tk.Listbox(
                motif_frame,
                selectmode=tk.MULTIPLE,
                yscrollcommand=scrollbar.set,
                height=20,
                bd=0,
                highlightthickness=0
            )

            self.motif_listbox.pack(fill="both",
                                    expand=True,
                                    padx=5,
                                    pady=5)
            scrollbar.config(command=self.motif_listbox.yview)

            self.all_motifs = [
                "ATG", "TATA", "CGCG", "AATT",
                "GGG", "TTT", "CCA", "GATA",
                "TGG", "CGT", "AAC", "GCGC",
                "TATATA", "ATAT"
            ]
            self.update_motif_list(self.all_motifs)



            # ===== PRAWA STRONA =====
            right_frame = ttk.Frame(main_frame)
            right_frame.pack(side="right", fill="both", expand=True)

            right_frame.rowconfigure(1, weight=1)
            right_frame.columnconfigure(0, weight=1)

            # ===== Panel sterowania =====
            control_frame = ttk.LabelFrame(right_frame, text="Dane wejściowe")
            control_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=10)

            control_frame.columnconfigure(1, weight=1)

            ttk.Label(control_frame, text="Plik:").grid(row=0, column=0, padx=5, pady=5)
            self.file_entry = ttk.Entry(control_frame)
            self.file_entry.grid(row=0, column=1, sticky="ew", padx=5)

            ttk.Button(control_frame,
                       text="Wybierz",
                       style="Primary.TButton",
                       command=self.browse_file).grid(row=0, column=2, padx=5)

            ttk.Label(control_frame, text="Accession:").grid(row=1, column=0, padx=5)
            self.accession_entry = ttk.Entry(control_frame)
            self.accession_entry.grid(row=1, column=1, sticky="ew", padx=5)

            ttk.Button(control_frame,
                       text="Analizuj",
                       style="Success.TButton",
                       command=self.analyze).grid(row=2, column=1, pady=10)

            ttk.Button(control_frame,
                       text="Eksport CSV",
                       style="Warning.TButton",
                       command=self.export_csv).grid(row=3, column=1, pady=5)

            ttk.Button(control_frame,
                       text="Eksport PDF",
                       style="Warning.TButton",
                       command=self.export_pdf).grid(row=4, column=1, pady=5)


            # =========================
            # ZAKŁADKI
            # =========================

            self.notebook = ttk.Notebook(right_frame)
            self.notebook.grid(row=1, column=0, sticky="nsew", padx=10, pady=10)

            # ===== Zakładka 1 – Tabela =====
            self.tab_table = ttk.Frame(self.notebook)
            self.notebook.add(self.tab_table, text="Tabela motywów")

            self.tab_table.rowconfigure(0, weight=1)
            self.tab_table.columnconfigure(0, weight=1)

            self.tree = ttk.Treeview(self.tab_table, show="headings")
            self.tree.grid(row=0, column=0, sticky="nsew")

            scrollbar_tree = ttk.Scrollbar(self.tab_table,
                                           orient="vertical",
                                           command=self.tree.yview)
            self.tree.configure(yscrollcommand=scrollbar_tree.set)
            scrollbar_tree.grid(row=0, column=1, sticky="ns")

            # ===== Zakładka 2 – Wykres =====
            self.tab_plot = ttk.Frame(self.notebook)
            self.notebook.add(self.tab_plot, text="Wykres motywów")

            self.tab_plot.rowconfigure(0, weight=1)
            self.tab_plot.columnconfigure(0, weight=1)

            self.plot_frame = ttk.Frame(self.tab_plot)
            self.plot_frame.grid(row=0, column=0, sticky="nsew")

            # ===== Zakładka 3 – Statystyki =====
            self.tab_stats = ttk.Frame(self.notebook)
            self.notebook.add(self.tab_stats, text="Statystyki")

            self.tab_stats.rowconfigure(0, weight=1)
            self.tab_stats.columnconfigure(0, weight=1)

            self.result_text = tk.Text(self.tab_stats)
            self.result_text.grid(row=0, column=0, sticky="nsew")

            # ===== Zakładka 4 – Heatmapa =====
            self.tab_heatmap = ttk.Frame(self.notebook)
            self.notebook.add(self.tab_heatmap, text="Heatmapa motywów w sekwencji DNA")

            self.tab_heatmap.rowconfigure(0, weight=1)
            self.tab_heatmap.columnconfigure(0, weight=1)

            self.heatmap_frame = ttk.Frame(self.tab_heatmap)
            self.heatmap_frame.grid(row=0, column=0, sticky="nsew")

            # ===== Zakładka CpG =====

            self.tab_cpg = ttk.Frame(self.notebook)
            self.notebook.add(self.tab_cpg, text="CpG Islands")

            self.tab_cpg.rowconfigure(0, weight=1)
            self.tab_cpg.columnconfigure(0, weight=1)

            self.cpg_tree = ttk.Treeview(self.tab_cpg, show="headings")

            self.cpg_tree["columns"] = ("Start", "End", "GC_content", "CpG_ratio")

            for col in self.cpg_tree["columns"]:
             self.cpg_tree.heading(col, text=col)

            self.cpg_tree.grid(row=0, column=0, sticky="nsew")

            scroll = ttk.Scrollbar(self.tab_cpg, orient="vertical",
                           command=self.cpg_tree.yview)

            self.cpg_tree.configure(yscrollcommand=scroll.set)
            scroll.grid(row=0, column=1, sticky="ns")

            # ===== Zakładka Mapa Genomu =====
            self.tab_genome = ttk.Frame(self.notebook)
            self.notebook.add(self.tab_genome, text="Mapa genomu")

            self.tab_genome.rowconfigure(0, weight=1)
            self.tab_genome.columnconfigure(0, weight=1)

            self.genome_frame = ttk.Frame(self.tab_genome)
            self.genome_frame.grid(row=0, column=0, sticky="nsew")



    # =============================
    def browse_file(self):
        filename = filedialog.askopenfilename()
        self.file_entry.delete(0, tk.END)
        self.file_entry.insert(0, filename)

    # =============================
    def analyze(self):

        selected_indices = self.motif_listbox.curselection()
        motifs = [self.motif_listbox.get(i) for i in selected_indices]

        if not motifs:
            messagebox.showerror("Błąd", "Wybierz przynajmniej jeden motyw!")
            return

        file_path = self.file_entry.get()
        accession = self.accession_entry.get()

        if file_path:
            self.sequence = load_sequence_from_file(file_path)
        elif accession:
            self.sequence = load_sequence_from_ncbi(accession)
            if self.sequence is None:
                messagebox.showerror("Błąd", "Błąd pobierania z NCBI")
                return
        else:
            messagebox.showerror("Błąd", "Wybierz źródło danych!")
            return

        self.results_df = segment_sequence_multiple(self.sequence, motifs)
        total_gc = gc_content(self.sequence)

        #sprawdzanie sekwencji
        if not validate_dna(self.sequence):
            from tkinter import messagebox

            messagebox.showerror(
                "Błąd sekwencji",
                "Plik zawiera znaki inne niż A, T, G, C.\nTo nie jest poprawna sekwencja DNA."
            )

            return


        #statystyki tekstowe
        self.result_text.delete("1.0", tk.END)
        self.result_text.insert(tk.END,
                                f"Długość sekwencji: {len(self.sequence)}\n")
        self.result_text.insert(tk.END,
                                f"GC-content (całość): {total_gc}%\n\n")
        self.notebook.select(self.tab_table)

        for motif in motifs:
            total = self.results_df[motif].sum()
            self.result_text.insert(tk.END,
                                    f"{motif}: {total} wystąpień\n")

        # tabeli
        self.tree["columns"] = ["Segment", "Start"] + motifs + ["GC_content_%"]

        for col in self.tree["columns"]:
            self.tree.heading(col, text=col)

        for row in self.tree.get_children():
            self.tree.delete(row)

        for _, row in self.results_df.iterrows():
            values = [row[col] for col in self.tree["columns"]]
            self.tree.insert("", "end", values=values)

        self.canvas = draw_plot(self.plot_frame,
                                self.results_df,
                                motifs,
                                self.canvas)
        #heatmapa
        draw_heatmap(self.heatmap_frame, self.results_df, motifs)

        # ===== CpG islands =====

        islands = find_cpg_islands(self.sequence)

        for row in self.cpg_tree.get_children():
            self.cpg_tree.delete(row)

        for island in islands:
            self.cpg_tree.insert("", "end", values=(
                island["Start"],
                island["End"],
                island["GC_content"],
                island["CpG_ratio"]
            ))

        draw_genome_map(
            self.genome_frame,
            self.results_df,
            motifs,
            islands
        )

    # =============================CSV===========================
    def export_csv(self):

        if self.results_df is None:
            messagebox.showerror("Błąd", "Najpierw wykonaj analizę!")
            return

        save_path = filedialog.asksaveasfilename(
            defaultextension=".csv"
        )

        if save_path:
            self.results_df.to_csv(save_path, index=False)
            messagebox.showinfo("Sukces", "Zapisano CSV")


    #======================================= PDF=================
    def export_pdf(self):

        from tkinter import filedialog
        from reportlab.lib.pagesizes import A4
        from reportlab.pdfgen import canvas

        if self.results_df is None:
            messagebox.showerror("Błąd", "Najpierw wykonaj analizę!")

            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF file", "*.pdf")]
        )

        if not file_path:
            return

        c = canvas.Canvas(file_path, pagesize=A4)

        y = 800

        c.setFont("Helvetica", 12)
        c.drawString(50, y, "DNA Motif Analysis Report")

        y -= 40
        c.drawString(50, y, f"Długość sekwencji: {len(self.sequence)}")

        y -= 30
        c.drawString(50, y, "Motywy:")

        y -= 20

        for motif in self.results_df.columns:
            if motif not in ["Segment", "Start", "GC_content_%"]:
                total = self.results_df[motif].sum()
                c.drawString(70, y, f"{motif}: {total}")
                y -= 20

        y -= 20
        c.drawString(50, y, "GC content segments")

        y -= 20

        for i, row in self.results_df.head(20).iterrows():

            text = f"Segment {row['Segment']}  GC% {row['GC_content_%']}"
            c.drawString(70, y, text)

            y -= 15

            if y < 50:
                c.showPage()
                y = 800

        c.save()