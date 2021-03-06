# SubPydence, A Multi-Layer, One-Dimensional Subsidence Model

## Introduction

SubPydence is a subsidence package built off of the modeling program, TTim, through which transient multi-layer flow with analytical elements can be modeled.
SubPydence consists of a library of Python scripts that add onto TTim, which is written in Python scripts with FORTRAN extensions for increased performance.

## TTim Installation

** Python versions:**

TTim requires **Python** 3.6 and can be installed by …

**Dependencies:**

TimML requires **Numpy** 1.12 (or higher) and **matplotlib** 2.0 (or higher).

**For base Python distributions:**

To install TTim, open a command prompt and type:

	pip install ttim

To update TTim type:

	pip install ttim —upgrade

To uninstall TTim type:

	pip uninstall ttim

**Testing installation:**

	python
	import ttim.ttimtest


## SubPydence Installation

**Python versions:**

SubPydence requires **Python** 3.6 

To install subpydence,

	git clone https://github.com/ecpease/subpydence.git

then in your script import sys and append a path to the subpydence directory

	import sys
	sys.path.append(os.path.join('somepath','subpydence'))
	import subpydence

**Dependencies:**

Additional dependencies require **ttim** version 0.4 (or higher), **os** version xxxx (or higher), and 
**sys** version xxx (or higher).

**For base Python distributions:**

To install SubPydence, open a command prompt and type:

	pip install subpydence

To update SubPydence type:

	pip install subpydence —upgrade

To uninstall SubPydence type:

	pip uninstall subpydence

**Testing installation:**

	python
	import ttim.subpydencetest

## Documentation

* The manual is available from the docs directory or can be viewed [here] (http://mbakker7.github.io/ttim/docs/builddocs/html/index.html).
* Example Notebooks are available from the notebooks directory on GitHub, or from [here] (https://github.commbakker7/ttim/tree/master/notebooks).

## Citation

We ask that if you found TTim (or SubPydence) useful to please include these citations in your papers:

* M. Bakker. 2013. Semi-analytic modeling of transient multi-layer flow with TTIm. Hydrogeology Journal, 21: 935-943.
* M. Bakker. 2013. Analytic modeling of transient multi-layer flow.  In: Advances in Hydrogeology, edited by P Mishra and K Kuhlman, Springer, Heidelberg, 95-114.
