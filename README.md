Computing Atomic Regulons
===============

Here we propose a new algorithm for computing ARs, which combines many of the advantages of the existing data-driven approaches, but integrates new evidence types including gene context and functional relationships to more quickly converge on a complete set of biologically meaningful ARs. Our algorithm is unique from other approaches in that it begins by constructing draft ARs using a combination of operon predictions and SEED subsystem technology.

Code to compute Atomic Regulons can be found here: https://github.com/jplfaria/atomic_regulons/blob/master/lib/Bio/KBase/atomic_regulons/atomic_regulonsImpl.pm


Test service:

1) Download RASTtk 1.3.0 available at: https://github.com/TheSEED/RASTtk-Distribution/releases/

The RASTtk/KBase environment is necessary for access to the SEED Subsystems and Functional roles

2) Launch RASTtk for to the RASTtk interactive shell

The RASTtk/KBase environment is necessary for access to the SEED Subsystems and Functional roles

3) Clone atomic_regulons repository:

git clone https://github.com/jplfaria/atomic_regulons.git

4) In the RASTtk interactive shell


