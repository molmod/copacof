# COPACOF: a computational protocol for the assignment of NMR resonances in covalent organic frameworks
<p align="center">
<img align="center" src="https://github.ugent.be/storage/user/1787/files/8b99b1c2-32bb-4f30-91d6-6dd24da78548">  
</p>

This repository contains data related to the following publication in ACS-JCTC: doi.org/10.1021/acs.jctc.3c01414. The folder `vasp/` contains input files for the structural optimization and NMR calculations of the covalent organic framework materials studied. The folder `spectrum/` contains the following Python scripts and files, which may be used to produce a computational NMR spectrum based on theoretical shielding data using the gs model as described in the main text:  
  (a) `classify_carbon.py`  
  (b) `gs_model.py`  
  (c) Input and output files as an example of the usage of the scripts.  

(a) **classify_carbon.py** Based on connectivity (using graphs), classify the carbon atoms in the COF into classes.  
    input files:  
      `system.extxyz`  --  periodic COF material  
      `node.xyz` -- TAPD node that will be searched for in the material  
      `linker.xyz` -- Me linker that will be searched for in the material  
    output files: txt files where the first column is the number of the carbon atom, which is a label derived from the index in the node/linker xyz. Next columns are the indices of the carbon atom in the COF material xyz that belong to that class. One class per row.  
      `linkers_classified.dat`  
      `nodes_classified.dat`  
      
 (b) **gs_model.py**  
    input files:  
      `shieldings` -- all shieldings in the COF system  
      `linkers_classified.dat`  
      `nodes_classified.dat`  
      `nodes_colordef.dat` -- color choice specific to example case. Change to whatever color codes desired.  
      `linkers_colordef.dat`  
      `spectrum.dat` -- the experimental spectrum to fit to  

Run: 
`python classify_carbon.py`  
`python gs_model.py mse`  

where mse is the choice for the residual function (others are implemented as well, see script). Additionally, boxplots (and other indicators) may be plotted for the data on the spectrum, which may be switched on in the script setting boxplot_on to value of 1.  

