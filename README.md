# teaching_code
Code that I have used for (or may be useful for) teaching.

Mostly targeted at simulations to explain genetics and evolution.

* `EvoSim.r`: written in Matlab originally in 2011 as one of the very first things I ever programmed, and rewrote it in R one year later for a one-day class I was teaching for the UW Science Education Partnership. Be gentle with it; I mostly didn't know what I was doing.
* `quant_gen_sim.R`: some prototype code for simulating phenotypes under different genetic architectures, with visualizations to demonstrate how different architectures affect population phenotypes. Needs some TLC but basically works for additive models. Handles different ploidies and numbers of contributing loci. Intended to illustrate the intuition behind Fisher (1918) resolving the competing Mendelist and Biometrician factions of genetics.
* `recombination.R`: a really brain-dead script that plots parental and recombinant chromosomes as an illustration.
* `exact_match.py`: answering a [question](https://bioinformatics.stackexchange.com/questions/14423/ncbi-blast-for-exact-match-of-a-short-sequence/14433#14433) on bioinformatics stackexchange site made me think that it would be useful to have a somewhat flexible sequence mapper based purely on exact matches (allowing Ns), using regex. 
