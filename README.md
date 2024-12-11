# Ambartsumyan function 

This repo is designed to those who might need to use Ambartsumyan function during their study or research in radiative transfer theory.
Examples of computed values are located in the folder tables, where files labeled with `_my_main` were generated using this code,
whereas other files are there for reference and comparison.

## How to use

All the necessary functions for computation are located in `ambcm_core.py` and you can use them simply by downloading this file
and writing `from ambcm_core import *` in your python script.

`ambcm_core` contains multiple versions of Ambartsumyan function catered for different use cases.
All of them share 2 arguments: z -- cosine of the reflection angle, $\lambda$ -- scattering albedo.

* `ambrcm(z, lamb)` -- main function for general use, significantly faster than the rest.
* `ambrcm_small(z, lamb)` -- special version for more precise values when `z` is low. Works slower.
* `ambrcm_mid(z, lamb)` -- special version for cases when z is big enough, yet still <= 1. Works slower.
* `ambrcm_auto(z, lamb)` -- automatically selects between `ambrcm_small` and `ambrcm_mid` depending on a value of `z`.

## Other python files

2 more files are added: `ambcm.py` allows use from console. `make_tables.py` creates reference tables for different values of 
z and $\lambda$. Tables take some time to be created, as `ambrcm_small` and `ambrcm_mid` are pretty slow due to having to compute multiple 
nestd integrals. 
All formulas and math are located in `calcfiambar.pdf` in Russian language.

