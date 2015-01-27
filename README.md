E. coli assembly
================================================================================

Pipeline
================================================================================

+ Connect reads using ABySS-mergepairs and Konnector
+ Assemble using ABySS-paireddbg
+ Scaffold using ABySS-scaffold and LINKS
+ Fill gaps using Sealer
+ Correct small-scale misassemblies using Scrubber

Sequencing data
================================================================================

Name: MiSeq v3: TruSeq Nano (E coli)
Read length: 301 bp
Fragment size: 600 bp

https://basespace.illumina.com/sample/3756762/Ecoli

Species: E. coli K-12 MG1655
Note: Not DH10B, which is the reference that Illumina used for alignment

See `/genesis/extscratch/btl/datasets/illuminaMiSeq_Fall2014/readme.txt`

> Ecoli_S1_L001_R1_001.fastq.gz
> Ecoli_S1_L001_R2_001.fastq.gz
>
> Downloaded from Illumina Base Space Oct 2014 (MiSeq data PE300 ~550bp fragment)
> ~2500X coverage of E.coli genome
>
> Split each file into 10, each expected to provide ~250-fold coverage of E.coli genome 

Relevant data and assemblies
================================================================================

+ `/genesis/extscratch/btl/paired-dbg/assembler-comparison/ecoli/data/merged`
+ `/genesis/extscratch/btl/paired-dbg/assembler-comparison/ecoli/abyss/full-assembly/k416`
+ `/genesis/extscratch/btl/paired-dbg/assembler-comparison/ecoli/abyss-paired-dbg/full-assembly/k32/d300`

File                              | Reads   | Proportion
----------------------------------|---------|--------
Ecoli_S1_L001_R1_001.fastq.5.gz   | 2239109 | 1
Ecoli_S1_L001_part5_merged.fastq  | 1065033 | 0.476
Ecoli_S1_L001_part5_reads_1.fastq | 1174076 | 0.524

Paired-end stats
================================================================================

```
Mapped 4396047 of 4478218 reads (98.2%)
Mapped 4373267 of 4478218 reads uniquely (97.7%)
Read 4478218 alignments
Mateless         0
Unaligned     7168  0.32%
Singleton    67835  3.03%
FR         2128682  95.1%
RF              96  0.00429%
FF              70  0.00313%
Different    35258  1.57%
Total      2239109
FR Stats mean: 589.4 median: 580 sd: 145.3 n: 2128441 min: 79 max: 1226 ignored: 337
               ____▁▁▂▂▃▄▅▆▇▇███████▇▆▆▅▅▄▃▃▂▂▂▁▁▁_____
```

Standard de Bruijn Graph
================================================================================

n    |n:500  |L50  |min  |N80    |N50     |N20     |E-size  |max     |sum      |name
---  |---    |---  |---  |---    |---     |---     |---     |---     |---      |---
109  |106    |14   |541  |57653  |106301  |209897  |129286  |313430  |4632995  |ecoli-unitigs.fa
66   |64     |10   |678  |90997  |179541  |257095  |172709  |358481  |4666371  |ecoli-contigs.fa
61   |59     |9    |678  |91400  |179541  |286168  |206848  |479753  |4666332  |ecoli-scaffolds.fa
78   |76     |13   |678  |80826  |113878  |212780  |143640  |314271  |4666332  |ecoli-scaftigs.fa

`/genesis/extscratch/btl/paired-dbg/assembler-comparison/ecoli/abyss/full-assembly/k416`

Paired de Bruijn Graph
================================================================================

Ben's assembly stats
------------------------------------------------------------

Uses se=merged-reads,unmerged-reads pe=unmerged-reads

n    |n:500  |L50  |min  |N80    |N50     |N20     |E-size  |max     |sum      |name
---  |---    |---  |---  |---    |---     |---     |---     |---     |---      |---
129  |122    |21   |594  |32940  |82237   |116108  |82922   |194094  |4607684  |ecoli-k32d300-unitigs.fa
89   |84     |16   |594  |59062  |107008  |176372  |112474  |222606  |4674904  |ecoli-k32d300-contigs.fa
69   |64     |10   |594  |80785  |179604  |240882  |168635  |359235  |4670672  |ecoli-k32d300-scaffolds.fa

`/genesis/extscratch/btl/paired-dbg/assembler-comparison/ecoli/abyss-paired-dbg/full-assembly/k32/d300`

Shaun's assembly stats
------------------------------------------------------------

Uses se=merged-reads,unmerged-reads pe=original-reads

n    |n:500  |L50  |min  |N80    |N50     |N20     |E-size  |max     |sum      |name
---  |---    |---  |---  |---    |---     |---     |---     |---     |---      |---
129  |122    |21   |594  |32940  |82237   |116108  |82922   |194094  |4607684  |ecoli-unitigs.fa
88   |83     |15   |594  |66816  |107087  |176296  |117910  |236527  |4650047  |ecoli-contigs.fa
64   |59     |10   |594  |89070  |176689  |265709  |177721  |363847  |4645835  |ecoli-scaffolds.fa
86   |73     |12   |717  |66816  |133800  |216544  |144896  |311673  |4643574  |ecoli-scaftigs.fa

`/projects/btl/sjackman/ecoli/k364-K32`

Scaffolding without contigging
================================================================================

n    |n:500  |L50  |min  |N80    |N50     |N20     |E-size  |max     |sum      |name
---  |---    |---  |---  |---    |---     |---     |---     |---     |---      |---
129  |122    |21   |594  |32940  |82237   |116108  |82922   |194094  |4607684  |k364-K32/ecoli-unitigs.fa
88   |83     |15   |594  |66816  |107087  |176296  |117910  |236527  |4650047  |k364-K32/ecoli-contigs.fa
64   |59     |10   |594  |89070  |176689  |265709  |177721  |363847  |4645835  |k364-K32/ecoli-scaffolds.fa
79   |72     |10   |594  |80734  |176986  |238758  |167066  |356466  |4605653  |k364-K32-scaff/ecoli-scaffolds.fa
86   |73     |12   |717  |66816  |133800  |216544  |144896  |311673  |4643574  |k364-K32/ecoli-scaftigs.fa
129  |113    |20   |699  |40995  |83220   |133172  |88279   |207675  |4603759  |k364-K32-scaff/ecoli-scaftigs.fa

`/projects/btl/sjackman/ecoli`

Sealer
================================================================================
