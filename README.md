# QS-MIPS-Numerics
> Numerically solves the governing equations in https://www.biorxiv.org/content/10.1101/2023.04.01.535124v1.abstract. 

[![NPM Version][npm-image]][npm-url]
[![Build Status][travis-image]][travis-url]
[![Downloads Stats][npm-downloads]][npm-url]

This code implements a Galerkin finite element method with quadratic Lagrange elements and BDF2 time stepping. The implementation uses the open-source library oomph-lib (see https://oomph-lib.github.io/oomph-lib/doc/html/index.html). 

![](header.png)

## Installation

After installing oomph-lib, place the folders 1D, 2D, 1Dsteady in the user_drivers folder and rerun the autogen.sh script from the oomph-lib installation. The code can then be compiled with a call to make. 

## Usage example



## CHANGELOG


## Meta

Distributed under the GNU LESSER GENERAL PUBLIC LICENSE. See ``LICENSE`` for more information.
