# CiEj_2016H66
This repositery contains files to simulate CiEj surfactants as described in this paper :

Simulating Bilayers of Nonionic Surfactants with the GROMOS-Compatible 2016H66 Force Field.
Senac *et al.* (2017), *Langmuir*, **33**, 10225–10238 ([DOI](http://doi.org/10.1021/acs.langmuir.7b01348))

The CiEj parameters are based on the GROMOS-compatible 2016H66 force field ([DOI](http://doi.org/10.1021/acs.jctc.6b00187)).

The repository is organized in different sub_directories named by GROMACS versions. Each sub-dir contains files to simulate CiEj using the corresponding GROMACS version.

## Sub-Directory `gromacs-4.0.7`
For each surfactant, you'll find 2 sub-directories:
- `small` which contains small bilayers of 90 surfactants simulated with the twin-range cutoff (group scheme) and the reaction-field (RF) correction.
- `large` which contains large bilayers of 1024 surfactants simulated with PME.

NB: The 2016H66 force field has been parameterized using RF, however it gave instabilities for large bilayers. Thus we only put PME params for large bilayers. We also did some tests on small bilayers using PME, and we got very similar results compared to RF simulations (see the Langmuir paper).

NB2: In the Langmuir paper, all simulations were simulated using GROMACS 4.0.7, which is the last version supporting RF. For GROMACS versions >= 4.5.*, RF is no longer supported
([see RF issue](https://redmine.gromacs.org/issues/1400)). We thus recommend to only use PME for any version >= 4.5.*.

## Sub-Directory `gromacs-5.0.7`
For each surfactant, you'll find one sub-directory:
- `small` which contains small bilayers of 90 surfactants simulated with PME. This has been tested by us and we get similar results as with version 4.0.7 using PME.

## Sub-directory `search_pore_ciej`
This sub-directory contains the python script (+ dependancies) to calculate pores within a CiEj bilayer. This code has been initially written by Caroline Senac ([initial repository](https://github.com/csenac/search_pore_ciej)).


All these files are under licence Creative Commons Attribution - Partage dans les Mêmes Conditions 3.0 France (CC BY-SA 3.0 FR).
