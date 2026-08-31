This folder is for setting up separate nextflow runs of Fungi_BFD for the same input samples of genomes, same RNASeq already downloaded, to run different versions (1.8.17, 1.9.0-beta.10) of funannotate with generally the same version of nextflow funannotate.nf (we will have to make changes to the nextflow to be able to point to differnet installed versions and challeges with module load vs container).

I want to evaluate
* performance improvement with the new custom rust implementation of dependencies (trinity,EVM,PASA)
* performance differences in changes to 1.8.17 -> 1.9.0 code improvements in funannotate
* performance differences in module / conda env vs container

The nextflow system is here git@github.com:stajichlab/Fungi_BFD.git but there is also stajichlab/nf_funannotate1 that is aiming to be stand alone nextflow funannotate run that has more generalization so you could try it.

I need
* Develop the input dataset of a random set of (N=100 or so, this should be a starting script that can be used to add more (and not remove any already selected in an existing samples.csv) genomes where never more than 1 from same species, should be a balance of taxonomic groups (ratio can be at least 10 from ascomycota making sure at least 5 from Peziziomycotina; 5 from basidiomycota with at least 2 from Agaricomycotina, 1 from Pucciniomycotina, and 1 from Ustilaginomycotina; 3 from mucoromycota, 1 from Chytridiomycota and 1 from Blastocladiomycota) . no genome size should be > 75Mb to keep this simple
* the random subset could be preferred from RefSeq which will have annotation available in NCBI download. These start with GCF_ typically. No more than 2 genomes from same Genus should chosen.
* we will use already masked genomes in previous analysis so you can symlink these into input_clean_genomes from ../Fungi_BFD_runs/input_clean_genomes

questions
* is it better to just try to compare everything built from containers anyways and not try deal with all the local custom nextflow changes it will require to run nextflow with module system again just for these tests? (we could check out old version of Fungi_BFD/nextflow or can we rely on nf_funannotate1)
* Provide summary on runtimes and memory usage based on nextflow or other profiling
* Provide summary on gene content prediction and assess differences in predictions (we could perhaps match against RefSeq GFF for genomes which have that)

Future
* also profile the impact of the pick representative tools in Fungi_BFD in ability to reuse prediction parameters to measure
  - speedup (no longer training augustus and genemark)
  - accuracy (perhaps better or same accuracy as individually trained on each strain)
