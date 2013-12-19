VQCDThermo
==========

A Mathematica package for computing thermodynamics (and more) in holographic Veneziano QCD. See [arXiv:1112.1261](http://arxiv.org/abs/1112.1261), [arXiv 1210.4516](http://arxiv.org/abs/1210.4516) and especially [arXiv 1312.5199](http://arxiv.org/abs/1312.5199) for background on VQCD and the type of computations this package does.

What this code does
-------------------

This code takes care of much of the boilerplate in doing numerical computations in the VQCD model at finite temperature and finite chemical potential. The lowest level work, that is, setting up the boundary conditions, solving the differential equations, fixing the various scales and so on, is contained in a fully automated form. The computation of the thermodynamic parameter scans is present in an incomplete form, and the thermodynamic analysis is currently not contained here. The intent is, however, to publish these in the near future, once the code used in [arxiv 1312.5199](http://arxiv.org/abs/1312.5199) has been brought to a form suitable for release.

Getting Started
===============

Installation
------------

The sub-directory VQCDThermo (not the similarly named top-level repository directory!) needs to be in Mathematica's search path. The usual place to install the package can be found as follows:

1.    Open a new Mathematica notebook. Type $UserBaseDirectory to an empty cell, and evaluate it. It should print out       a directory on your local filesystem.
2.    Navigate to that directory on you computer. It should have a sub-directory Applications.
3.    Copy the sub-directory VQCDThermo from the repository into the Applications sub-directory

The file should now be discoverable by Mathematica. Verify this by evaluating Needs["VQCDThermo`"] in a notebook cell. It should not produce any output, especially you should not see a Get:noopen warning, if everything went right.

You are then ready to move to the next step.

Computing black hole background solutions and extracting basic thermo observables
---------------------------------------------------------------------------------

First of all, read [arXiv 1312.5199](http://arxiv.org/abs/1312.5199), and make sure you have a basic understanding on the purpose of the code, and the various observables. Then open Examples/BasicTutorial.nb, which shows how to carry out the basic numerical computations.

Computing thermodynamics
------------------------

In order to compute full thermodynamics in this model, you need to scan over the boundary parameters \lambda_h and \tilde{n}. This is handled, including automatic parallelization, file saving, finding to correct region in the parameter space (for some potentials), bookkeeping and so on, by the functions in VQCDThermo.m. However, at the moment, this code is still somewhat untested (there are some known problems), and the corresponding analysis stage is not yet in shape to be published.

The computations in [arXiv 1312.5199](http://arxiv.org/abs/1312.5199) were done using a less automated (i.e. requiring manual parameter setup and manual parallelization) and more complicated predecessor to this code (the parts in VQCDCore.m are pretty much unchanged), and a somewhat laborious semi-manual analysis stage. Once we finish the automatic code, it will be released here.

Getting help
------------

You can always contact me by e-mail, and I'll do my best to help you out with any problems concerning the code. If your problem is clearly a bug in the code, you can also just leave an issue in https://github.com/timoalho/VQCDThermo/issues.

TODO:
=====

- [ ]   Finish the automated thermo and analysis code
- [ ]   Write better documentation
- [ ]   Include initial conditions in the core solvers to also construct the vacuum solutions
- [ ]   Set up a method in the core solvers to allow this code to be used for IHQCD calculations also (they are a            special case of calculations in VQCD anyway). 
