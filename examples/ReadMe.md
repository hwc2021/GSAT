
# Test Dataset Description

This directory contains a test dataset for the complete mitochondrial genome assembly pipeline of rice sample **NH002**.

## File Structure

- **`run.sh`**: Shell script for the entire workflow.

> **Note**: Whole-genome sequencing data files are too large to include here. For the full test dataset (including ONT reads and Illumina PE FASTQ files), please visit Zenodo:  
> https://doi.org/10.5281/zenodo.18730781

## Per-Step Configuration and Results

To save storage, only essential result files are kept; intermediate files have been removed.

### 1. graphShort Pipeline

- Config file: `short.conf`
- Output directory: `NH002_short/`

#### Additional step for removing non-organelle sequences

- **clean graph** → `NH002_short/og.filtered.clean.gfa`

### 2. graphLong Pipeline

- Config file: `long.conf`
- Output directory: `NH012_long/`

#### Additional step for merging potential linear paths

- **merged graph** → `NH012_long/mrg.filtered.merge.gfa`

### 3. graphSimplification Pipeline

- Config file: `sim.conf`
- Output directory: `NH012_sim/`

### 4. graphCorrection Pipeline

- Config file: `corr.conf`
- Output directory: `NH012_corr/`

### 5. Final Merge

- Final MMG: `NH012_corr/mmg.corrCtg.merge.gfa`

