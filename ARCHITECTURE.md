# GeneWeaver – System Architecture

## 1. Overview

GeneWeaver follows a modular architecture. Each module handles a separate responsibility, including DNA alignment, GPU processing, distributed candidate screening, FASTA parsing, PAM detection, and dashboard interaction.

## 2. Architecture Diagram

```text
                 USER
                   |
                   v
          TEXTUAL DASHBOARD
                   |
          +--------+--------+
          |        |        |
          v        v        v
        FASTA    ALIGNMENT  BENCHMARK
        INPUT     ENGINE    ENGINE
          |        |
          v        v
        PARSER   +----------------+
                 |                |
                 v                v
             CPU ENGINE       GPU ENGINE
                 |                |
                 +-------+--------+
                         |
                         v
                  ALIGNMENT RESULTS
                         |
             +-----------+-----------+
             |           |           |
             v           v           v
          DASK        PAM         CSV EXPORT
        SCREENING   DETECTION
             |
             v
       RANKED CANDIDATES
```

## 3. Core Modules

### cpu_alignment.py

Implements Needleman–Wunsch global alignment.

Responsibilities:
- Compare two DNA sequences.
- Identify matches and mismatches.
- Calculate insertions and deletions.
- Calculate similarity percentage.
- Generate alignment scores.

### gpu_alignment.py

Provides CUDA-based GPU processing.

Responsibilities:
- Execute supported sequence comparisons on the GPU.
- Process alignment workloads using CUDA.
- Support different sequence lengths through GPU matrix processing.
- Compare GPU results against CPU results.

Note: The equal-length GPU comparison path uses position-wise comparison, whereas the CPU engine uses global alignment. These methods are not interchangeable for every sequence pair.

### distributed_alignment.py

Uses Dask for distributed task scheduling.

Responsibilities:
- Schedule candidate alignment tasks.
- Process multiple candidates independently.
- Rank candidates by similarity and alignment score.
- Return structured screening results.

The current parallelism is across candidate sequences rather than across individual cells of one dynamic-programming alignment matrix.

### parser.py

Handles FASTA file processing.

Responsibilities:
- Read FASTA files.
- Extract sequence identifiers.
- Extract DNA sequences.
- Validate and normalize sequence data.

### pam_scoring.py

Implements prototype PAM detection.

Responsibilities:
- Detect SpCas9 NGG PAM motifs.
- Search forward and reverse orientations.
- Calculate prototype PAM scores.
- Return detected PAM locations.

PAM detection and scoring are computational prototypes, not validated predictions of biological editing efficiency.

### dashboard/app.py

Provides the interactive Textual dashboard.

Responsibilities:
- Load FASTA files.
- Select target and candidate sequences.
- Run CPU and GPU comparisons.
- Perform Dask candidate screening.
- Display benchmark information.
- Display alignment results.
- Export screening results to CSV.

## 4. Processing Workflow

1. User launches the dashboard.
2. User loads a FASTA file.
3. FASTA parser extracts DNA sequences.
4. User selects a target and candidate.
5. CPU and GPU alignment engines process supported comparisons.
6. Dask schedules candidate screening tasks.
7. Candidate results are ranked.
8. PAM detection identifies matching motifs.
9. Results are displayed in the dashboard.
10. Candidate screening results can be exported to CSV.

## 5. Testing Architecture

Pytest validates the project's modules.

The test suite covers:
- CPU alignment.
- GPU alignment.
- Distributed screening.
- FASTA parsing.
- PAM detection.
- Result validation.

Latest recorded result: **23 passed, 2 warnings**.

## 6. Future Improvements

- Advanced GPU optimization.
- More efficient alignment algorithms.
- Improved dashboard visualizations.
- Larger FASTA dataset support.
- More comprehensive benchmarking.
- Additional CRISPR PAM models.
- Integration with validated biological reference datasets.

---

**Project:** GeneWeaver  
**Domain:** Bioinformatics, GPU Computing, Distributed Computing  
**Language:** Python