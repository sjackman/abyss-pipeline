---
title: ABySS assembly of S. cerevisiae
author: Shaun Jackman
---

Data
================================================================================

[S288C Reference Genome](http://www.yeastgenome.org/download-data/sequence)

Oxford Nanopore
------------------------------------------------------------

+ [Oxford Nanopore Sequencing and de novo Assembly of a Eukaryotic Genome](http://schatzlab.cshl.edu/data/nanocorr/)
+ Oxford Nanopore raw reads: [W303_ONT_Raw_reads.fa.gz](http://labshare.cshl.edu/shares/schatzlab/www-data/nanocorr/W303_ONT_Raw_reads.fa.gz)
+ Nanocorr corrected reads: [W303_ONT_Nanocorr_Corrected_reads.fa.gz](http://labshare.cshl.edu/shares/schatzlab/www-data/nanocorr/W303_ONT_Nanocorr_Corrected_reads.fa.gz)
+ Nanocorr polished assembly: [W303_ONT_Assembly_CA_polished.fa.gz](http://labshare.cshl.edu/shares/schatzlab/www-data/nanocorr/W303_ONT_Assembly_CA_polished.fa.gz)

Longest Repeat
================================================================================

```
repeat-match -n500 R64.fa
Genome Length = 230218   Used 302109 internal nodes
Long Exact Matches:
   Start1     Start2    Length
    26322     205049r      817
```

```
exonerate --bestn 2 --showalignment false R64.fa R64.fa
vulgar: X 8 18847 + IX 25 18864 + 94159 M 18839 18839
vulgar: IX 25 18864 + X 8 18847 + 94159 M 18839 18839
vulgar: VIII 532690 543613 + I 212088 223011 + 54462 M 10923 10923
vulgar: I 212088 223011 + VIII 532690 543613 + 54462 M 10923 10923
vulgar: IV 7721 17328 + X 738078 728471 - 47603 M 9607 9607
vulgar: XII 460554 468943 + XII 451417 459806 + 41900 M 8389 8389
vulgar: XIV 21 6616 + VII 1090940 1084345 - 32957 M 6595 6595
vulgar: VII 1084345 1090940 + XIV 6616 21 - 32957 M 6595 6595
vulgar: XV 1084588 1091142 + XII 1071623 1078177 + 32635 M 6554 6554
vulgar: XVI 0 6528 + VII 1090870 1084342 - 32640 M 6528 6528
vulgar: III 84807 90786 + VI 137909 143888 + 29652 M 5979 5979
vulgar: II 29642 35604 + IV 519648 513686 - 29144 M 5962 5962
vulgar: VI 137912 143858 + VII 811446 817392 + 29712 M 5946 5946
vulgar: V 443396 449320 + XVI 62375 56451 - 29521 M 5924 5924
vulgar: XIII 196331 202242 + V 449320 443409 - 29456 M 5911 5911
vulgar: XI 661401 663351 + III 8366 6416 - 9444 M 1950 1950
vulgar: Mito 21213 21603 + Mito 41523 41913 + 1347 M 390 390
```
