## Como instalar recount
## https://bioconductor.org/packages/recount
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("recount")

## Checar nuestra versión de BioC
BiocManager::version()

## Checar que tengamos una instalación
## válida de Bioconductor
BiocManager::valid()
