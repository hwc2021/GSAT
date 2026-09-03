source /public/home/heweichuang123/miniforge3/bin/activate software

#fasterq-dump -3 -e 16 -p SRR13453724.sra

#seqkit sample -s 11 -p 0.1 -o SRR13453724.s01.1.fq.gz SRR13453724_1.fastq
#seqkit sample -s 11 -p 0.1 -o SRR13453724.s01.2.fq.gz SRR13453724_2.fastq

#seqkit sample -s 11 -p 0.25 ../NH012.ONT.mt_mapped.fasta > ../NH012.ONT.mt_mapped.s025.fasta

# Short reads
read2_1=SRR13453724.s01.1.fq.gz
read2_2=SRR13453724.s01.2.fq.gz

# Long reads
read3=NH012.mt.s025.fasta

# graphShort pipeline
[ -d NH012_short ] && rm -r NH012_short
gsat graphShort -conf short.conf -cpu 12

# Clean the OG file with a new script "anchor_subgraph.pl". You can also do this with Bandage.
perl ../scripts/anchor_subgraph.pl -g NH012_short/og.filtered.gfa -r mh63.mmg.m0.fasta -o NH012_short/og.filtered.clean.gfa

# graphLong pipeline
[ -d NH012_long ] && rm -r NH012_long
gsat graphLong -conf long.conf
# It is recommended to conduct a merging operation after each graphLong pipeline.
gsat graphMerge -g NH012_long/mrg.filtered.gfa -o NH012_long/mrg.filtered

# Conduct the graphMap for the subsequent graphClip process. The graphMap–graphClip workflow will be implemented as a pipeline in the next update.
gsat graphMap -a -d -g NH012_long/mrg.filtered.merge.gfa -r $read3 -maxOffset1 150 -maxCombDis 150 -maxEdgeSize1 150 -maxEdgeSize2 150 -minIden 0.8 -minCovofRead 0.8 -minCovbyPath 0.8 -o NH012_mrg -minimap2 ont

# graphClip
gsat graphClip -g NH012_long/mrg.filtered.merge.gfa -mtReadList NH012_mrg.mt.reads -o NH012_mrg_clip5 -minPathNum 5 -b7File NH012_mrg.part2.reads.b7 -pafFile NH012_mrg.part1.reads.paf

# graphSimplification pipeline
[ -d NH012_sim ] && rm -r NH012_sim
gsat graphSimplification -conf sim.conf

# graphCorrection pipeline
[ -d NH012_corr ] && rm -r NH012_corr
gsat graphCorrection -conf corr.conf

# graphMerge 
gsat graphMerge -g NH012_corr/mmg.corrCtg.gfa -o NH012_corr/mmg.corrCtg

# Note: We recommend performing an additional check for each MMG using the graphSimplification pipeline to ensure that all links are well supported by long reads and that all repeat-mediated recombinant structures have been thoroughly evaluated. If not, additional graphSimplification–graphCorrection–graphMerge steps can be performed to resolve these issues.
