## Pasos para crear este RStudio proyect
## con este archivo de notas
usethis::create_project("~/rnaseq_2024_notas")
usethis::use_r("01-notas.R")

## Ejercicio 2
usethis::use_r("02-visualizar-mtcars.R")

## Ligar este proyecto a GitHub
usethis::use_git()
usethis::use_github()

## Notas para intro a BioC
usethis::use_r("03-intro-Bioconductor.R")

## Notas sobre SummarizedExperiment
usethis::use_r("04_SummarizedExperiment.R")

## Notas sobre recount3
usethis::use_r("05_recount3.R")
