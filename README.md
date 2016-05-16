meta-sweeper
============

By expressing microbial community composition and the details of metagenomic
experiments in parametric form, meta-sweeper aims to permit the assessment of
analytical methods under variation in both community and experiment.

Standard (local) execution
```bash
nextflow -C sweep.config run sweep.nf
```

SGE execution
```bash
nextflow -C sweep.config run sweep.nf --profile sge
```

```bash
nextflow -C sweep.comfig run sweep.nf --profile pbspro
```

[Darling Lab](http://darlinglab.org/)

[i3 institute, UTS
Sydney](http://www.uts.edu.au/research-and-teaching/our-research/ithree-institute)
