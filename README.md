ONTotalRecaller
===============
Authors: Nicola Prezza, Bud Mishra, Giuseppe Narzisi
mail: nicolapr@gmail.com

## Overview

this software is still under development!


ONTotalRecaller (ONTRC) is a tool designed to align Oxford Nanopore events on a reference genome. ONTRC performs both base-call of ON raw signals and alignment at the same time, being able to (given a good reference) correct signal errors that are hard to detect using standard base callers (e.g. homopolymers)

## Download

Since ONTotalRecaller includes extern git repositories as submodules, clone it using the --recursive option:

> git clone --recursive https://github.com/nicolaprezza/ONTotalRecaller

## Install extern libraries

ONTotalRecaller requires the libraries sdsl and hdf5 to be installed in your system. Read below for more info on how to do this.

### HDF5

Note: we suggest to download HDF5 version 1.8 (at the moment, version 1.9 creates problems in the parsing of FAST5 files). Download one of the .tar.bz2 files at ftp://ftp.hdfgroup.uiuc.edu/pub/outgoing/hdf5/snapshots/v18/ , unpack it, enter the unpacked directory and execute:

> ./configure --prefix=/usr/local/hdf5

> make

> sudo make install

For more informations, see https://www.hdfgroup.org/ftp/HDF5/current/src/unpacked/release_docs/INSTALL

### SDSL

Clone the sdsl library with

> git clone https://github.com/simongog/sdsl-lite

Enter the directory sdsl-lite and execute

> ./install.sh

sdsl headers and libraries will be automatically installed in ~/include and ~/lib

## Compile

At the moment, only debug code is present. Wait until the tool is finshed:)


ONTRC uses cmake to automatically generate a Makefile. Create a build directory, run cmake, and compile ONTRC as follows (from ONTotalRecaller directory):

> mkdir build

> cd build

> cmake ..

> make
