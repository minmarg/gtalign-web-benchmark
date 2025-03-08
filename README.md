# GTalign-web benchmark scripts

GTalign-web benchmark scripts for benchmarking GTalign-web, a web 
[implementation](https://github.com/minmarg/gtalign-web-backend) of 
[GTalign](https://github.com/minmarg/gtalign_alpha), 
against other protein structure alignment webservers.

## Description

The `commands-test.sh` provides commented commands for running and evaluating all tested webservers. 

The `queries` directory contains 100 query structures used to benchmark webservers.

The `get_tmscores_dali_server.pl`, `get_tmscores_foldseek_server.pl`, and `get_tmscores_gtalign_web.pl` 
scripts evaluate alignments produced by the DALI and Foldseek servers and GTalign-web, respectively, 
by calculating TM-scores and GDT\_TS scores.

The `TMscore_from_alignment` program required by these scripts to calculate GDT\_TS scores is 
available at [GTalign-evaluation](https://github.com/minmarg/gtalign-evaluation). 
These scripts also use `getchain.py` to extract chains from PDB files.

The `plot_TMscores_GDTTSs.sh` script generates the main becnhmark results figure.

These directories contain the raw data required for the plotting script to generate graphs:

  * `Dali-server-processed.gdtts.out`
  * `Dali-server-processed.out`
  * `Foldseek-server-3diaa-processed.gdtts.out`
  * `Foldseek-server-3diaa-processed.out`
  * `Foldseek-server-tmalign-processed.gdtts.out`
  * `Foldseek-server-tmalign-processed.out`
  * `gtalign_web_speed13_prescore04-processed.gdtts.out`
  * `gtalign_web_speed9_prescore04-processed.gdtts.out`

GTalign-web's benchmark results are available as three completed jobs using a speed setting of 13 (speed=13):

  * [Part1](https://bioinformatics.lt/comer/gtalign/results/benchmark13_p1)
  * [Part2](https://bioinformatics.lt/comer/gtalign/results/benchmark13_p2)
  * [Part3](https://bioinformatics.lt/comer/gtalign/results/benchmark13_p3)

and using a speed setting of 9 (speed=9):

  * [Part1](https://bioinformatics.lt/comer/gtalign/results/benchmark9_p1)
  * [Part2](https://bioinformatics.lt/comer/gtalign/results/benchmark9_p2)
  * [Part3](https://bioinformatics.lt/comer/gtalign/results/benchmark9_p3)

The results of GTalign-web's case studies from UniRef30 searches are available at the following link:

  * [Casestudies](https://bioinformatics.lt/comer/gtalign/results/example_uniref30)


