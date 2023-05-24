# GSAT - Graph-based Sequence Assembly Toolkit

An effective graph-based toolkit to assembly plant organelle genomes into simple and accurate master graphs, consisting of many graph-based tools for processing the genome assembly results and high throughout sequencing data.

## Installation

Just download the file and uncompress it to where you want to install in. The executable file `gsat` could be found in the `bin` directory.

    #install by using git
    git clone https://github.com/hwc2021/GSAT.git
    cd GSAT/bin
    chmod a+x gsat
    
    #install by downloading the source codes
    #put the source code file "GSAT-main.zip" where you want to install in
    unzip GSAT-main.zip
    cd GSAT-main/bin
    chmod a+x gsat
    
Edit the `.bashrc` file to add this directory to your environment variable `PATH`. For example, in a typical Linux system:

    vi ~/.bashrc
    #add the next line to the end of .bashrc file ("#" should be removed when paste the next line to the file)
    #export PATH=$PATH:/your/path/GSAT/bin
    source ~/.bashrc

## Requirement

This toolkit required several necessary tools which should be installed in the `PATH`:

- [**blastn**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [**minimap2**](https://github.com/lh3/minimap2)
- [**mafft**](https://mafft.cbrc.jp/alignment/software/)
- [**SPAdes**](https://github.com/ablab/spades) [used only in graphShort pipeline]

## Usage

    gsat <command> [options]

    Commands:
    -- Functions
       graphFilt            filter the assembly graph with different params
       graphMap             conduct graph mapping to detect mapped paths in a graph for query sequence
       graphCorr            correct the sequences in a graph by using long reads. HIFI reads is recommanded.
       graphSimplify        simplify the graph based on supported mapped paths of long reads.
       rmOverlap            remove the overlaping regions from a graph

    -- Pipelines
       graphShort           generate a Organelle Graph from a raw graph of de novo assembly
       graphLong            generate a Mitochondrial Rough Graph from a OG
       graphSimplification  generate a Mitochondrial Rough Master Graph from a MRG
       graphCorrection      generate a Mitochondrial Master Graph from a MRMG

    -- Information
       help                 print a brief help information
       man                  print a complete help document
       version              print the version information
       
Detailed usage information for each command could be found by using the command without any option.

## Examples for using pipelines

All pipeline commands can be easily applied with a specified configure file, which provides all necessary params to the pipelines of GSAT. The detailed params for each pipeline can be found in the example configure file `example.conf`. Notably, each pipeline reads only valid options from the configure file, i.e., params for different pipelines can be either put together in a single configure file or put in different files. 

- graphShort

        gsat graphShort -conf example.conf

This will generate a raw graph of de novo genome assembly based on Illumina paired-end reads, and then produce a Organelle Graph (OG). The OG will be saved as `og.filtered.gfa` in the output dir. Please note that the params for graphShort should be adjusted for different species, so the resulted OG should be checked and edited in [Bandage](https://github.com/rrwick/Bandage) before it can be used for the next pipeline.

- graphLong

        gsat graphLong -conf example.conf

This will generate a Mitochondrial Rough Graph (MRG) from a OG based on graph-mapping between OG and long reads. The MRG will be saved as `mrg.filtered.gfa` in the output dir.

- graphSimplification

        gsat graphSimplification -conf example.conf

This will generate a Mitochondrial Rough Master Graph (MRMG) from a MRG. The MRMG will be saved as `mrmg.simplified.gfa` in the output dir.

- graphCorrection

        gsat graphCorrection -conf example.conf

This will generate a Mitochondrial Master Graph (MMG) from a MRMG. The MMG will be saved as `mmg.corrCtg.gfa` in the output dir.

## Version

GSAT version 1.10 (2022-09-30)

## Contact

  Wenchuang He
  nongke2@163.com / hewenchuang@caas.cn
  AGIS, CAAS

## Citation

He, W., Xiang, K., Chen, C., Wang, J., and Wu, Z. (2023). Master graph: an essential integrated assembly model for the plant mitogenome based on a graph-based framework. Brief Bioinform 2410.1093/bib/bbac522.

## Feedback
Any comments, bug reports, and suggestions are very welcomed. They will help us to further improve GSAT. If you have any troubles running GSAT, please leave your comments and bug reports at our GitHub [**issues page**](https://github.com/hwc2021/GSAT/issues) or sent it via e-mail.
