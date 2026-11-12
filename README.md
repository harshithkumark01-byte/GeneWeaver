# 🧬 GeneWeaver
### GPU-Accelerated CRISPR DNA Alignment and Candidate Screening Engine

GeneWeaver is a Python-based DNA sequence alignment and candidate screening system designed to explore CPU processing, GPU acceleration, and distributed computing using Dask.

The project provides DNA sequence comparison, FASTA file parsing, candidate ranking, PAM detection, benchmarking, and CSV result export through an interactive terminal dashboard.

## 🚀 Key Features

- **CPU Alignment:** Needleman–Wunsch global DNA sequence alignment.
- **GPU Acceleration:** CUDA-based DNA sequence comparison using NVIDIA GPUs.
- **Distributed Computing:** Dask-based parallel candidate screening.
- **FASTA Processing:** Load and parse DNA sequences from FASTA files.
- **Candidate Ranking:** Rank candidate sequences using alignment similarity and score.
- **PAM Detection:** Identify SpCas9 NGG PAM sequences.
- **Benchmarking:** Compare CPU and GPU execution times.
- **CSV Export:** Export candidate screening results.
- **Interactive Dashboard:** Terminal-based interface built using Textual and Rich.

## 🛠️ Technologies Used

| Technology | Purpose |
|---|---|
| Python | Core programming language |
| Numba CUDA | GPU acceleration |
| NVIDIA CUDA | GPU computing |
| Dask | Distributed task scheduling |
| Textual | Interactive terminal dashboard |
| Rich | Terminal formatting |
| NumPy | Numerical computing |
| Pytest | Automated testing |

## 📁 Project Structure

```text
GeneWeaver/
│
├── src/
│   └── geneweaver/
│       ├── cpu_alignment.py
│       ├── gpu_alignment.py
│       ├── distributed_alignment.py
│       ├── parser.py
│       ├── pam_scoring.py
│       └── dashboard/
│           └── app.py
│
├── data/
│   ├── sample_genome.fasta
│   └── alignment_pair.fasta
│
├── tests/
│
├── results/
│   └── candidate_screening_results.csv
│
├── README.md
└── ARCHITECTURE.md
```

## ⚙️ Installation

Clone the repository:

```bash
git clone YOUR_GITHUB_REPOSITORY_URL
```

Navigate into the project:

```bash
cd GeneWeaver
```

Create a virtual environment:

```bash
python -m venv .venv
```

Activate it on Windows:

```powershell
.venv\Scripts\Activate.ps1
```

Install dependencies:

```bash
pip install -r requirements.txt
```

For GPU acceleration, an NVIDIA GPU and a compatible CUDA environment are required.

## ▶️ Running GeneWeaver

Set the Python source path:

```powershell
$env:PYTHONPATH="src"
```

Launch the dashboard:

```powershell
python -m geneweaver.dashboard.app
```

## 🧪 Running Tests

Execute the complete test suite:

```powershell
python -m pytest tests/ -v
```

Latest test result:

**23 passed, 2 warnings**

## 📊 Performance Benchmarking

GeneWeaver compares CPU and GPU execution times across different DNA sequence sizes.

Example benchmark results:

| Sequence Size | CPU Time | GPU Time | Speedup |
|---|---:|---:|---:|
| 10,000 | 0.78 ms | 2.61 ms | 0.30x |
| 100,000 | 10.62 ms | 3.72 ms | 2.85x |
| 1,000,000 | 78.06 ms | 19.24 ms | 4.06x |
| 5,000,000 | 427.25 ms | 94.09 ms | 4.54x |

*These are example measurements from the project's previous benchmark run. Actual performance depends on hardware, workload, and implementation.*

## 🧬 Scientific Scope

GeneWeaver is an educational and computational prototype for DNA sequence alignment and candidate screening.

Its alignment scores and PAM-based rankings are computational indicators. They do not establish experimental CRISPR editing efficiency, biological safety, or clinical suitability.

## 👨‍💻 Project Goals

- Explore GPU-accelerated sequence processing.
- Implement global DNA sequence alignment.
- Apply distributed computing to candidate screening.
- Compare CPU and GPU performance.
- Develop an interactive bioinformatics tool.

---

**GeneWeaver – Exploring DNA Alignment Through Parallel Computing.**