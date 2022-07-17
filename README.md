# GSAT - Graph-based Sequence Assembly Toolkit

An effective graph-based toolkit to assembly plant organelle genomes into simple and accurate master graphs, consisting of many graph-based tools for processing the genome assembly results and high throughout sequencing data.

## Installation

Just download the file and uncompress it to where you want to install in. The executable file `gsat` could be found in the `bin` directory.

## Requirement

This toolkit required several necessary tools which should be installed in the `PATH`:

- [**blastn**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [**minimap2**](https://github.com/lh3/minimap2)
- [**mafft**](https://mafft.cbrc.jp/alignment/software/) 

## Usage

    ./bin/gsat <command> [options]

    Commands:
    -- Functions
       graphFilt            filter the assembly graph with different params
       graphMap             conduct graph mapping to detect mapped paths in a graph for query sequence
       graphCorr            correct the sequences in a graph by using long reads. HIFI reads is recommanded.
       graphSimplify        simplify the graph based on supported mapped paths of long reads.
       rmOverlap            remove the overlaping regions from a graph

    -- Pipelines (incoming)
       graphShort           generate a Organelle Graph from a raw graph of de novo assembly
       graphLong            generate a Mitochondrial Rough Graph from a OG
       graphSimplification  generate a Mitochondrial Rough Master Graph from a MRG
       graphCorrection      generate a Mitochondrial Master Graph from a MRMG

    -- Information
       help                 print a brief help information
       man                  print a complete help document
       version              print the version information

## Version

GSAT version 1.00 (2022-07-17)

## Author

  Wenchuang He
  nongke2@163.com
  AGIS, CAAS
