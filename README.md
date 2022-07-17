# GSAT - Graph-based Sequence Assembly Toolkit

An effective toolkit to assembly plant organelle genomes into simple and accurate master graphs, consisting of many graph-based tools for processing the genome assembly results and high throughout sequencing data.

## SYNOPSIS

    gsat <command> [options]

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

## VERSION

GSAT version 1.00 (2022-07-17)

## REQUIREMENTS

This toolkit required several necessary tools which could be found in the PATH variable:

    blastn
    minimap2
    mafft 

## AUTHOR

  Wenchuang He
  nongke2@163.com
  AGIS, CAAS
