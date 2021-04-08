# Human Brain Metabolic Network

**Scope:** This repository contains all the necessary code to generate the centrality and optimality analysis of the neuron-astrocyte metabolic network. The network was derived from Lewis

## Proyect structure

### Directories

- **source** includes main source code for this proyect. Eg. `.py`, `.r`, `.rmd`. Includes also the `.sbatch` and `.sh` files used for cluster excecution in SLURM. 
- **data** includes `.json` COBRA models. 
- **results** includes network models `.gephx`, output tables `.csv` and `.xsls`. 
- **doc** includes the source code for the paper. It's mostly `.rmd` files rendered with RMarkdown. 
    - **img** includes generated figures `.png`, `.svg`, etc. 

### Branches

- `main` contains the main development branch, reused 
- `brain-astrocyte` contains development and specific code related to the 

### Flow del proyecto

![](./doc/img/flow_codigo.drawio.svg)

### .gitignore

Excludes the `tmp.*` files, and `tmp/` folder, which contains large datasets available in AWS S3. 
ALso excludes RStudio files (`.RData`. `.Rhistory`, ...), and `.pem` key files. 

-----

## References



