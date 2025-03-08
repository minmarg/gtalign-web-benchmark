# GTalign-web benchmark scripts

GTalign-web benchmark scripts for benchmarking GTalign-web, a web 
[implementation](https://github.com/minmarg/gtalign-web-backend) of 
[GTalign](https://github.com/minmarg/gtalign_alpha), 
against other protein structure alignment webservers.

## Description

The `commands-test.sh` provides commented commands for running and evaluating all tested webservers. 

The `get_tmscores_dali_server.pl`, `get_tmscores_foldseek_server.pl`, and `get_tmscores_gtalign_web.pl` 
scripts evaluate alignments produced by the DALI and Foldseek servers and GTalign-web, respectively, 
by calculating TM-scores and GDT\_TS scores.

The `TMscore_from_alignment` program required by these scripts to calculate GDT\_TS scores is 
available at [GTalign-evaluation](https://github.com/minmarg/gtalign-evaluation). 
These scripts also use `getchain.py` to extract chains from PDB files.

The `plot_TMscores_GDTTSs.sh` script generates the main becnhmark results figure.

The `queries` directory contains 100 query structures used to benchmark webservers.

