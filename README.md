Computing Atomic Regulons
===============

Here we propose a new algorithm for computing ARs, which combines many of the advantages of the existing data-driven approaches, but integrates new evidence types including gene context and functional relationships to more quickly converge on a complete set of biologically meaningful ARs. Our algorithm is unique from other approaches in that it begins by constructing draft ARs using a combination of operon predictions and SEED subsystem technology.

Code to compute Atomic Regulons can be found here: https://github.com/jplfaria/atomic_regulons/blob/master/lib/Bio/KBase/atomic_regulons/atomic_regulonsImpl.pm


Test service
-------------

1) Download RASTtk 1.3.0 available at: https://github.com/TheSEED/RASTtk-Distribution/releases/

Note: 
The RASTtk/KBase environment is necessary for access to the SEED Subsystems and Functional roles

2) Launch RASTtk distribution to prompt the RASTtk interactive shell

3) Clone atomic_regulons repository:

```
git clone https://github.com/jplfaria/atomic_regulons.git

```

4) In the RASTtk interactive shell run the following cmd in atomic_regulons/test-service : 

e.g., for Escherichia coli data, run:

```
perl -I ../lib testARserviceImpl.pl "kb|g.0" e.coli_expression.tab

```
Parameters:
- Genome ID: "kb|g.0"
- Expression Data: e.coli_expression.tab (provided in /test-service)

Note:
Search for genome ID for genome of interest here: https://narrative.kbase.us/functional-site/#/search/?q=ecoli%20k12&category=genomes

5) Output Atomic regulons are available at:

/test-service/ar.out