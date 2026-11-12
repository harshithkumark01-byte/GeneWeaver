import asyncio
import time
import csv
import subprocess

from pathlib import Path

from textual.app import App, ComposeResult
from textual.widgets import (
    Header,
    Footer,
    Static,
    ProgressBar,
    Button,
    Input,
    Select
)
from textual.containers import Vertical, Horizontal, VerticalScroll

from rich.text import Text

from geneweaver.cpu_alignment import compare_sequences
from geneweaver.gpu_alignment import compare_sequences_gpu

from geneweaver.distributed_alignment import (
    distributed_alignment,
    distributed_screen
)

from geneweaver.pam_scoring import score_pam_sites
from geneweaver.parser import parse_fasta


PROJECT_ROOT = Path(__file__).resolve().parents[3]


class GeneWeaverDashboard(App):

    TITLE = "GeneWeaver - CRISPR Alignment Engine"

    CSS = """
    Screen {
        background: #101827;
    }

    VerticalScroll {
        height: 1fr;
        scrollbar-size: 1 1;
    }

    .panel {
        border: round #00bcd4;
        padding: 1;
        margin: 1;
        height: auto;
    }

    .panel-title {
        color: #00e5ff;
        text-style: bold;
        margin-bottom: 1;
    }

    Static {
        color: white;
        height: auto;
    }

    Input, Select {
        margin: 1;
        border: round #00bcd4;
    }

    #button_bar {
        height: 7;
        width: 100%;
        align: center middle;
        background: #101827;
    }

    Button {
        width: 25;
        margin: 1;
    }

    #run_button {
        background: #00897b;
        color: white;
    }

    #refresh_button {
        background: #3949ab;
        color: white;
    }

    #load_button {
        background: #0277bd;
        color: white;
    }

    #screen_button {
        background: #8e24aa;
        color: white;
    }

    #export_button {
        background: #2e7d32;
        color: white;
    }

    ProgressBar {
        margin: 1;
    }
    """

    def __init__(self):

        super().__init__()

        self.target_sequence = None
        self.candidate_sequence = None

        self.target_id = None
        self.candidate_id = None

        self.fasta_loaded = False
        self.fasta_records = {}

        # CSV EXPORT DATA

        self.last_screening_result = None
        self.last_screening_target = None
        self.last_screening_time = None
        self.last_screening_report = ""

    # ------------------------------------------------
    # DASHBOARD LAYOUT
    # ------------------------------------------------

    def compose(self) -> ComposeResult:

        yield Header()

        with VerticalScroll():

            yield Static(
                "GENEWEAVER GPU-ACCELERATED CRISPR ALIGNMENT ENGINE",
                classes="panel-title"
            )

            # FASTA INPUT

            with Vertical(classes="panel"):

                yield Static(
                    "FASTA SEQUENCE INPUT",
                    classes="panel-title"
                )

                yield Static(
                    "Load a FASTA file containing multiple DNA sequences."
                )

                yield Input(
                    value="data/sample_genome.fasta",
                    placeholder="Enter FASTA file path",
                    id="fasta_path"
                )

                yield Static(
                    "FASTA STATUS",
                    id="fasta_status"
                )

                yield Static(
                    "TARGET SEQUENCE",
                    classes="panel-title"
                )

                yield Select(
                    [],
                    prompt="Select target sequence",
                    allow_blank=True,
                    id="target_select"
                )

                yield Static(
                    "CANDIDATE SEQUENCE",
                    classes="panel-title"
                )

                yield Select(
                    [],
                    prompt="Select candidate sequence",
                    allow_blank=True,
                    id="candidate_select"
                )

                yield Static(
                    "Load a FASTA file to see available sequences.",
                    id="sequence_information"
                )

            # GPU INFORMATION

            with Vertical(classes="panel"):

                yield Static(
                    "GPU INFORMATION",
                    classes="panel-title"
                )

                yield Static(
                    "Loading GPU information...",
                    id="gpu_info"
                )

            # ALIGNMENT RESULTS

            with Vertical(classes="panel"):

                yield Static(
                    "GLOBAL ALIGNMENT RESULTS",
                    classes="panel-title"
                )

                yield Static(
                    "Waiting for alignment...",
                    id="alignment_results"
                )

            # ALIGNED SEQUENCES

            with Vertical(classes="panel"):

                yield Static(
                    "NEEDLEMAN-WUNSCH DNA ALIGNMENT",
                    classes="panel-title"
                )

                yield Static(
                    "Aligned sequences will appear here.",
                    id="aligned_sequences"
                )

            # MISMATCH DETAILS

            with Vertical(classes="panel"):

                yield Static(
                    "MISMATCH AND INDEL ANALYSIS",
                    classes="panel-title"
                )

                yield Static(
                    "Waiting...",
                    id="mismatch_total"
                )

                yield Static(
                    "Run alignment to see mismatch details.",
                    id="mismatch_details"
                )

            # COLORED DNA

            with Vertical(classes="panel"):

                yield Static(
                    "COLORED DNA SEQUENCE COMPARISON",
                    classes="panel-title"
                )

                yield Static(
                    "Green = Match | Red = Mismatch | Yellow = Gap"
                )

                yield Static(
                    "TARGET DNA:",
                    id="target_dna"
                )

                yield Static(
                    "CANDIDATE DNA:",
                    id="candidate_dna"
                )

                yield Static(
                    "COMPARISON:",
                    id="dna_comparison"
                )

            # PROGRESS

            with Vertical(classes="panel"):

                yield Static(
                    "ALIGNMENT PROGRESS",
                    classes="panel-title"
                )

                yield ProgressBar(
                    total=100,
                    show_eta=False,
                    id="progress"
                )

                yield Static(
                    "Waiting for alignment...",
                    id="progress_text"
                )

            # PERFORMANCE

            with Vertical(classes="panel"):

                yield Static(
                    "CPU VS GPU PERFORMANCE",
                    classes="panel-title"
                )

                yield Static(
                    "Waiting...",
                    id="performance"
                )

            # DASK

            with Vertical(classes="panel"):

                yield Static(
                    "DASK DISTRIBUTED RESULTS",
                    classes="panel-title"
                )

                yield Static(
                    "Waiting...",
                    id="dask_results"
                )

            # DASK PARALLEL SCREENING

            with Vertical(classes="panel"):

                yield Static(
                    "DASK PARALLEL DNA CANDIDATE SCREENING",
                    classes="panel-title"
                )

                yield Static(
                    "Screen all candidate sequences against the selected target."
                )

                yield Static(
                    "Waiting for candidate screening...",
                    id="screening_results"
                )

            # VALIDATION

            with Vertical(classes="panel"):

                yield Static(
                    "CPU VS GPU VS DASK VALIDATION",
                    classes="panel-title"
                )

                yield Static(
                    "Waiting for validation...",
                    id="validation_results"
                )

            # PAM

            with Vertical(classes="panel"):

                yield Static(
                    "CRISPR PAM ANALYSIS",
                    classes="panel-title"
                )

                yield Static(
                    "Waiting...",
                    id="pam_results"
                )

            # BENCHMARK

            with Vertical(classes="panel"):

                yield Static(
                    "BENCHMARK VISUALIZATION",
                    classes="panel-title"
                )

                yield Static(
                    "Loading benchmark...",
                    id="benchmark_chart"
                )

        # BUTTON BAR

        with Horizontal(id="button_bar"):

            yield Button(
                "LOAD FASTA",
                id="load_button"
            )

            yield Button(
                "RUN ALIGNMENT",
                id="run_button"
            )

            yield Button(
                "SCREEN CANDIDATES",
                id="screen_button"
            )

            yield Button(
                "EXPORT CSV",
                id="export_button",
                disabled=True
            )

            yield Button(
                "REFRESH BENCHMARK",
                id="refresh_button"
            )

        yield Footer()

    # ------------------------------------------------
    # GPU INFORMATION
    # ------------------------------------------------

    def get_gpu_information(self):

        try:

            result = subprocess.run(
                [
                    "nvidia-smi",
                    "--query-gpu=name,memory.used,memory.total,utilization.gpu",
                    "--format=csv,noheader,nounits"
                ],
                capture_output=True,
                text=True,
                timeout=5
            )

            if result.returncode == 0:

                values = result.stdout.strip().split(",")

                return (
                    f"GPU: {values[0].strip()}\n"
                    f"Memory: {values[1].strip()} / "
                    f"{values[2].strip()} MB\n"
                    f"GPU Utilization: {values[3].strip()}%"
                )

        except Exception:
            pass

        return "GPU information unavailable."

    # ------------------------------------------------
    # FASTA LOADING
    # ------------------------------------------------

    def load_fasta_file(self):

        file_path = self.query_one(
            "#fasta_path",
            Input
        ).value.strip()

        if not file_path:

            raise ValueError("Please enter a FASTA file path.")

        path = Path(file_path)

        if not path.is_absolute():

            path = PROJECT_ROOT / path

        path = path.resolve()

        records = parse_fasta(path)

        if len(records) < 2:

            raise ValueError(
                "FASTA file must contain at least two sequences."
            )

        self.fasta_records = records

        options = [
            (sequence_id, sequence_id)
            for sequence_id in records.keys()
        ]

        target_select = self.query_one(
            "#target_select",
            Select
        )

        candidate_select = self.query_one(
            "#candidate_select",
            Select
        )

        target_select.set_options(options)
        candidate_select.set_options(options)

        sequence_ids = list(records.keys())

        target_select.value = sequence_ids[0]
        candidate_select.value = sequence_ids[1]

        self.fasta_loaded = True

        self.update_selected_sequences()

        self.query_one(
            "#fasta_status",
            Static
        ).update(
            f"FASTA LOADED SUCCESSFULLY!\n"
            f"File: {path.name}\n"
            f"Total sequences: {len(records)}"
        )

        self.query_one(
            "#progress_text",
            Static
        ).update(
            "FASTA loaded. Select sequences and run alignment."
        )

    # ------------------------------------------------
    # UPDATE SELECTED SEQUENCES
    # ------------------------------------------------

    def update_selected_sequences(self):

        target_id = self.query_one(
            "#target_select",
            Select
        ).value

        candidate_id = self.query_one(
            "#candidate_select",
            Select
        ).value

        if not isinstance(target_id, str):
            return

        if not isinstance(candidate_id, str):
            return

        if target_id not in self.fasta_records:
            return

        if candidate_id not in self.fasta_records:
            return

        self.target_id = target_id
        self.candidate_id = candidate_id

        self.target_sequence = self.fasta_records[target_id]
        self.candidate_sequence = self.fasta_records[candidate_id]

        target_length = len(self.target_sequence)
        candidate_length = len(self.candidate_sequence)

        self.query_one(
            "#sequence_information",
            Static
        ).update(
            f"Available FASTA sequences: {len(self.fasta_records)}\n\n"
            f"Target: {self.target_id}\n"
            f"Target Length: {target_length:,} bases\n\n"
            f"Candidate: {self.candidate_id}\n"
            f"Candidate Length: {candidate_length:,} bases\n\n"
            f"Length Difference: {abs(target_length - candidate_length):,} bases\n\n"
            "Status: READY FOR GLOBAL ALIGNMENT\n"
            "Different-length sequences are supported."
        )

    # ------------------------------------------------
    # DROPDOWN EVENTS
    # ------------------------------------------------

    def on_select_changed(self, event: Select.Changed):

        if event.select.id in (
            "target_select",
            "candidate_select"
        ):

            if self.fasta_loaded:
                self.update_selected_sequences()

    # ------------------------------------------------
    # ALIGNED DNA VISUALIZATION
    # ------------------------------------------------

    def display_dna_comparison(
        self,
        aligned_target,
        aligned_candidate
    ):

        target_text = Text("TARGET: ")
        candidate_text = Text("CANDIDATE: ")
        comparison_text = Text("COMPARE: ")

        display_limit = 150

        target_window = aligned_target[:display_limit]
        candidate_window = aligned_candidate[:display_limit]

        for base1, base2 in zip(
            target_window,
            candidate_window
        ):

            if base1 == base2:

                target_text.append(base1, style="bold green")
                candidate_text.append(base2, style="bold green")
                comparison_text.append("|", style="green")

            elif base1 == "-":

                target_text.append("-", style="bold yellow")
                candidate_text.append(base2, style="bold yellow")
                comparison_text.append("I", style="yellow")

            elif base2 == "-":

                target_text.append(base1, style="bold yellow")
                candidate_text.append("-", style="bold yellow")
                comparison_text.append("D", style="yellow")

            else:

                target_text.append(base1, style="bold red")
                candidate_text.append(base2, style="bold red")
                comparison_text.append("X", style="bold red")

        self.query_one("#target_dna", Static).update(target_text)
        self.query_one("#candidate_dna", Static).update(candidate_text)
        self.query_one("#dna_comparison", Static).update(comparison_text)

    # ------------------------------------------------
    # BENCHMARK CHARTS
    # ------------------------------------------------

    def load_benchmark_charts(self):

        csv_file = PROJECT_ROOT / "benchmark_results.csv"

        if not csv_file.exists():

            self.query_one(
                "#benchmark_chart",
                Static
            ).update(
                "Benchmark file not found.\n"
                "Run: python -m geneweaver.benchmark"
            )

            return

        try:

            with csv_file.open(
                "r",
                newline="",
                encoding="utf-8"
            ) as file:

                rows = list(csv.DictReader(file))

            if not rows:
                raise ValueError("Benchmark CSV is empty.")

            chart = "CPU EXECUTION TIME (ms)\n"

            max_cpu = max(
                float(row["CPU_Time_ms"])
                for row in rows
            )

            for row in rows:

                value = float(row["CPU_Time_ms"])
                bar = "█" * max(1, int(value / max_cpu * 24))

                chart += (
                    f"{int(row['DNA_Size']):,}: "
                    f"{bar} {value:.2f} ms\n"
                )

            chart += "\nGPU EXECUTION TIME (ms)\n"

            max_gpu = max(
                float(row["GPU_Time_ms"])
                for row in rows
            )

            for row in rows:

                value = float(row["GPU_Time_ms"])
                bar = "█" * max(1, int(value / max_gpu * 24))

                chart += (
                    f"{int(row['DNA_Size']):,}: "
                    f"{bar} {value:.2f} ms\n"
                )

            chart += "\nGPU SPEEDUP\n"

            max_speedup = max(
                float(row["GPU_Speedup"])
                for row in rows
            )

            for row in rows:

                value = float(row["GPU_Speedup"])
                bar = "█" * max(1, int(value / max_speedup * 24))

                chart += (
                    f"{int(row['DNA_Size']):,}: "
                    f"{bar} {value:.2f}x\n"
                )

            self.query_one(
                "#benchmark_chart",
                Static
            ).update(chart)

        except Exception as error:

            self.query_one(
                "#benchmark_chart",
                Static
            ).update(
                f"Benchmark error: {error}"
            )

    # ------------------------------------------------
    # STARTUP
    # ------------------------------------------------

    def on_mount(self):

        self.query_one(
            "#gpu_info",
            Static
        ).update(
            self.get_gpu_information()
        )

        self.load_benchmark_charts()

        try:

            self.load_fasta_file()

            self.notify("Multi-sequence FASTA loaded!")

        except Exception as error:

            self.query_one(
                "#fasta_status",
                Static
            ).update(
                f"FASTA loading failed: {error}"
            )

    # ------------------------------------------------
    # BUTTON EVENTS
    # ------------------------------------------------

    async def on_button_pressed(self, event: Button.Pressed):

        if event.button.id == "load_button":

            try:

                self.load_fasta_file()

                self.last_screening_result = None
                self.last_screening_target = None
                self.last_screening_time = None
                self.last_screening_report = ""

                self.query_one(
                    "#export_button",
                    Button
                ).disabled = True

                self.notify("FASTA file loaded successfully!")

            except Exception as error:

                self.fasta_loaded = False
                self.fasta_records = {}

                self.query_one(
                    "#fasta_status",
                    Static
                ).update(
                    f"FASTA ERROR: {error}"
                )

                self.notify(
                    str(error),
                    severity="error",
                    timeout=8
                )

        elif event.button.id == "run_button":

            await self.run_alignment()

        elif event.button.id == "screen_button":

            await self.run_candidate_screening()

        elif event.button.id == "export_button":

            self.export_screening_csv()

        elif event.button.id == "refresh_button":

            self.load_benchmark_charts()
            self.notify("Benchmark refreshed!")

    # ------------------------------------------------
    # MAIN ALIGNMENT
    # ------------------------------------------------

    async def run_alignment(self):

        button = self.query_one("#run_button", Button)
        progress = self.query_one("#progress", ProgressBar)

        button.disabled = True

        try:

            if not self.fasta_loaded:
                self.load_fasta_file()

            self.update_selected_sequences()

            target = self.target_sequence
            candidate = self.candidate_sequence

            if target is None or candidate is None:
                raise ValueError("Please select both DNA sequences.")

            progress.update(progress=0)

            self.query_one(
                "#progress_text",
                Static
            ).update("Running CPU Needleman-Wunsch alignment...")

            cpu_start = time.perf_counter()

            cpu_result = await asyncio.to_thread(
                compare_sequences,
                target,
                candidate
            )

            cpu_time = (time.perf_counter() - cpu_start) * 1000

            progress.update(progress=25)

            self.query_one(
                "#progress_text",
                Static
            ).update("Initializing NVIDIA GPU...")

            await asyncio.to_thread(
                compare_sequences_gpu,
                "ATGC" * 1000,
                "ATGC" * 1000
            )

            self.query_one(
                "#progress_text",
                Static
            ).update("Running GPU global alignment...")

            gpu_start = time.perf_counter()

            gpu_result = await asyncio.to_thread(
                compare_sequences_gpu,
                target,
                candidate
            )

            gpu_time = (time.perf_counter() - gpu_start) * 1000

            progress.update(progress=50)

            self.query_one(
                "#progress_text",
                Static
            ).update("Running Dask global alignment...")

            dask_start = time.perf_counter()

            dask_result = await asyncio.to_thread(
                distributed_alignment,
                target,
                candidate,
                10000
            )

            dask_time = (time.perf_counter() - dask_start) * 1000

            progress.update(progress=70)

            self.query_one(
                "#progress_text",
                Static
            ).update("Analyzing CRISPR PAM sites...")

            pam_results = await asyncio.to_thread(
                score_pam_sites,
                target
            )

            progress.update(progress=85)

            def same_alignment(first, second):

                return (
                    first["aligned_target"] == second["aligned_target"]
                    and first["aligned_candidate"] == second["aligned_candidate"]
                    and first["alignment_score"] == second["alignment_score"]
                )

            cpu_gpu_match = same_alignment(cpu_result, gpu_result)
            cpu_dask_match = same_alignment(cpu_result, dask_result)

            validation = cpu_gpu_match and cpu_dask_match

            self.query_one(
                "#alignment_results",
                Static
            ).update(
                f"Target: {self.target_id}\n"
                f"Candidate: {self.candidate_id}\n\n"
                f"Original Target Length: {len(target):,}\n"
                f"Original Candidate Length: {len(candidate):,}\n"
                f"Aligned Length: {cpu_result['alignment_length']:,}\n\n"
                f"Matches: {cpu_result['matches']:,}\n"
                f"Mismatches: {len(cpu_result['mismatches']):,}\n"
                f"Insertions: {len(cpu_result['insertions']):,}\n"
                f"Deletions: {len(cpu_result['deletions']):,}\n"
                f"Similarity: {cpu_result['similarity']}%\n"
                f"Alignment Score: {cpu_result['alignment_score']}"
            )

            self.query_one(
                "#aligned_sequences",
                Static
            ).update(
                f"CPU ALIGNED TARGET:\n"
                f"{cpu_result['aligned_target']}\n\n"
                f"CPU ALIGNED CANDIDATE:\n"
                f"{cpu_result['aligned_candidate']}\n\n"
                f"GPU Alignment Score: {gpu_result['alignment_score']}\n"
                f"Dask Alignment Score: {dask_result['alignment_score']}"
            )

            mismatch_text = "MISMATCHES:\n"

            if cpu_result["mismatches"]:

                for item in cpu_result["mismatches"][:20]:

                    mismatch_text += (
                        f"Position: {item['position']} | "
                        f"Target: {item['target']} | "
                        f"Candidate: {item['candidate']}\n"
                    )

            else:
                mismatch_text += "None\n"

            mismatch_text += "\nINSERTIONS:\n"

            if cpu_result["insertions"]:

                for item in cpu_result["insertions"][:20]:

                    mismatch_text += (
                        f"Position: {item['position']} | "
                        f"Inserted Base: {item['candidate']}\n"
                    )

            else:
                mismatch_text += "None\n"

            mismatch_text += "\nDELETIONS:\n"

            if cpu_result["deletions"]:

                for item in cpu_result["deletions"][:20]:

                    mismatch_text += (
                        f"Position: {item['position']} | "
                        f"Deleted Base: {item['target']}\n"
                    )

            else:
                mismatch_text += "None\n"

            self.query_one(
                "#mismatch_total",
                Static
            ).update(
                f"Substitutions: {len(cpu_result['mismatches'])}\n"
                f"Insertions: {len(cpu_result['insertions'])}\n"
                f"Deletions: {len(cpu_result['deletions'])}\n"
                f"Total Differences: {cpu_result['mismatch_count']}"
            )

            self.query_one(
                "#mismatch_details",
                Static
            ).update(mismatch_text)

            self.display_dna_comparison(
                cpu_result["aligned_target"],
                cpu_result["aligned_candidate"]
            )

            progress.update(progress=90)

            speedup = cpu_time / gpu_time if gpu_time > 0 else 0

            self.query_one(
                "#performance",
                Static
            ).update(
                f"CPU Global Alignment: {cpu_time:.4f} ms\n"
                f"GPU Execution: {gpu_time:.4f} ms\n"
                f"GPU Speedup: {speedup:.2f}x\n\n"
                "Note: GPU time includes execution overhead."
            )

            self.query_one(
                "#dask_results",
                Static
            ).update(
                f"Dask Total Time: {dask_time:.4f} ms\n"
                f"Dask Matches: {dask_result['matches']:,}\n"
                f"Dask Mismatches + Gaps: {dask_result['mismatch_count']:,}\n"
                f"Dask Similarity: {dask_result['similarity']}%\n"
                f"Dask Alignment Score: {dask_result['alignment_score']}\n"
                f"Dask Tasks: {dask_result['chunks_processed']}\n"
                f"Execution Mode: {dask_result['execution_mode']}"
            )

            self.query_one(
                "#validation_results",
                Static
            ).update(
                f"CPU vs GPU: {'MATCH' if cpu_gpu_match else 'DIFFERENT'}\n"
                f"CPU vs Dask: {'MATCH' if cpu_dask_match else 'DIFFERENT'}\n\n"
                f"Overall Validation: {'SUCCESS' if validation else 'WARNING'}\n\n"
                "Validation compares aligned sequences and alignment scores."
            )

            forward_count = sum(
                site["strand"] == "+"
                for site in pam_results
            )

            reverse_count = sum(
                site["strand"] == "-"
                for site in pam_results
            )

            pam_text = (
                f"Total PAM Sites: {len(pam_results):,}\n"
                f"Forward Strand: {forward_count:,}\n"
                f"Reverse Strand: {reverse_count:,}\n\n"
                "PAM Details:\n"
            )

            for site in pam_results[:10]:

                pam_text += (
                    f"Position: {site['position']} | "
                    f"PAM: {site['pam']} | "
                    f"Strand: {site['strand']} | "
                    f"Score: {site['score']}\n"
                )

            pam_text += (
                "\nNote: PAM scores are prototype scores, "
                "not experimentally validated predictions."
            )

            self.query_one(
                "#pam_results",
                Static
            ).update(pam_text)

            progress.update(progress=100)

            status = (
                "SUCCESS: CPU, GPU and Dask alignments match!"
                if validation
                else "WARNING: Alignment results differ. Check validation."
            )

            self.query_one(
                "#progress_text",
                Static
            ).update(f"Alignment completed! {status}")

            if validation:
                self.notify("CPU, GPU and Dask validation successful!")
            else:
                self.notify(
                    "Alignment results differ. Check validation.",
                    severity="warning",
                    timeout=8
                )

        except Exception as error:

            self.query_one(
                "#progress_text",
                Static
            ).update(
                f"Alignment failed: {type(error).__name__}: {error}"
            )

            self.notify(
                f"Alignment failed: {error}",
                severity="error",
                timeout=8
            )

        finally:
            button.disabled = False

    # ------------------------------------------------
    # DASK PARALLEL CANDIDATE SCREENING
    # ------------------------------------------------

    async def run_candidate_screening(self):

        button = self.query_one("#screen_button", Button)
        export_button = self.query_one("#export_button", Button)

        output = self.query_one("#screening_results", Static)

        button.disabled = True
        export_button.disabled = True

        self.last_screening_result = None
        self.last_screening_target = None
        self.last_screening_time = None
        self.last_screening_report = ""

        try:

            if not self.fasta_loaded:
                self.load_fasta_file()

            self.update_selected_sequences()

            target_id = self.target_id
            target = self.target_sequence

            if not target:
                raise ValueError("Please select a target sequence.")

            candidates = [
                {
                    "id": sequence_id,
                    "sequence": sequence
                }
                for sequence_id, sequence in self.fasta_records.items()
                if sequence_id != target_id
            ]

            if not candidates:
                raise ValueError("No candidate sequences available.")

            output.update(
                "Running Dask parallel candidate screening...\n"
                f"Target: {target_id}\n"
                f"Candidates: {len(candidates)}"
            )

            start_time = time.perf_counter()

            screening_result = await asyncio.to_thread(
                distributed_screen,
                target,
                candidates
            )

            elapsed_time = (
                time.perf_counter() - start_time
            ) * 1000

            self.last_screening_result = screening_result
            self.last_screening_target = target_id
            self.last_screening_time = elapsed_time

            report = (
                "GENEWEAVER PARALLEL SCREENING RESULTS\n"
                "====================================\n\n"
                f"Target: {target_id}\n"
                f"Target Length: {len(target):,} bases\n"
                f"Candidates Processed: "
                f"{screening_result['candidates_processed']}\n"
                f"Scheduler: {screening_result['scheduler']}\n"
                f"Execution Mode: {screening_result['execution_mode']}\n"
                f"Total Screening Time: {elapsed_time:.4f} ms\n\n"
                "RANKED CANDIDATES\n"
                "-----------------\n\n"
            )

            for result in screening_result["results"]:

                report += (
                    f"Rank: {result['rank']}\n"
                    f"Candidate: {result['candidate_id']}\n"
                    f"Length: {result['candidate_length']:,} bases\n"
                    f"Matches: {result['matches']:,}\n"
                    f"Mismatches + Gaps: {result['mismatch_count']:,}\n"
                    f"Similarity: {result['similarity']}%\n"
                    f"Alignment Score: {result['alignment_score']}\n"
                    f"{'-' * 35}\n"
                )

            report += (
                "\nSCREENING COMPLETED SUCCESSFULLY!\n\n"
                "Note: Candidates are ranked using alignment "
                "similarity and alignment score. These rankings "
                "are computational comparisons, not experimentally "
                "validated CRISPR guide efficacy predictions."
            )

            # STORE THE REPORT FOR CSV EXPORT DISPLAY

            self.last_screening_report = report

            output.update(report)

            export_button.disabled = False

            self.notify("Dask candidate screening completed!")

        except Exception as error:

            output.update(
                f"SCREENING FAILED:\n"
                f"{type(error).__name__}: {error}"
            )

            self.notify(
                f"Screening failed: {error}",
                severity="error",
                timeout=8
            )

        finally:
            button.disabled = False

    # ------------------------------------------------
    # EXPORT SCREENING RESULTS TO CSV
    # ------------------------------------------------

    def export_screening_csv(self):

        if not self.last_screening_result:

            self.notify(
                "Run candidate screening first!",
                severity="warning"
            )
            return

        try:

            results = self.last_screening_result["results"]

            if not results:
                raise ValueError("No screening results available.")

            output_directory = PROJECT_ROOT / "results"

            output_directory.mkdir(
                parents=True,
                exist_ok=True
            )

            csv_path = (
                output_directory /
                "candidate_screening_results.csv"
            )

            columns = [
                "Target_ID",
                "Target_Length",
                "Rank",
                "Candidate_ID",
                "Candidate_Length",
                "Matches",
                "Substitutions",
                "Insertions",
                "Deletions",
                "Total_Differences",
                "Similarity_Percent",
                "Alignment_Score",
                "Scheduler",
                "Execution_Mode",
                "Screening_Time_ms"
            ]

            with csv_path.open(
                "w",
                newline="",
                encoding="utf-8"
            ) as file:

                writer = csv.DictWriter(
                    file,
                    fieldnames=columns
                )

                writer.writeheader()

                for result in results:

                    writer.writerow({

                        "Target_ID": self.last_screening_target,

                        "Target_Length": (
                            self.last_screening_result["target_length"]
                        ),

                        "Rank": result["rank"],

                        "Candidate_ID": result["candidate_id"],

                        "Candidate_Length": result["candidate_length"],

                        "Matches": result["matches"],

                        "Substitutions": len(result["mismatches"]),

                        "Insertions": len(result["insertions"]),

                        "Deletions": len(result["deletions"]),

                        "Total_Differences": result["mismatch_count"],

                        "Similarity_Percent": result["similarity"],

                        "Alignment_Score": result["alignment_score"],

                        "Scheduler": self.last_screening_result["scheduler"],

                        "Execution_Mode": self.last_screening_result["execution_mode"],

                        "Screening_Time_ms": self.last_screening_time
                    })

            output = self.query_one(
                "#screening_results",
                Static
            )

            output.update(
                self.last_screening_report
                + "\n\nCSV EXPORTED SUCCESSFULLY!\n"
                + f"File: {csv_path}"
            )

            self.notify(
                "CSV exported successfully!",
                timeout=8
            )

        except Exception as error:

            self.notify(
                f"CSV export failed: {error}",
                severity="error",
                timeout=8
            )


if __name__ == "__main__":

    GeneWeaverDashboard().run()