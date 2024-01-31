## Código de https://lcolladotor.github.io/rnaseq_LCG-UNAM_2024/objetos-de-bioconductor-para-datos-de-expresi%C3%B3n.html

## Lets build our first SummarizedExperiment object
library("SummarizedExperiment")
## ?SummarizedExperiment

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6
## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
## Información de nuestros genes
rowRanges <- GenomicRanges::GRanges(
    rep(c("chr1", "chr2"), c(50, 150)),
    IRanges::IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
    strand = sample(
        x = c("+", "-"),
        size = 200,
        replace = TRUE
    ),
    feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))
## Información de nuestras muestras
colData <- DataFrame(
    Treatment = rep(c("ChIP", "Input"), 3),
    row.names = LETTERS[1:6]
)
## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
    assays = SimpleList(counts = counts),
    rowRanges = rowRanges,
    colData = colData
)

## Exploremos el objeto resultante
rse

## Algunos ejemplos de comandos para explorar nuestro
## objeto de RangedSummarizedExperiment (RSE)
dim(rse)
assayNames(rse)
head(assay(rse))
rowRanges(rse)
mcols(rowRanges(rse))
rowData(rse)
rowData(rse)$feature_id[3]
colData(rse)
colData(rse)$Treatment[4]
rse$Treatment[4]


## Comando 1
rse[1:2, ]
## Comando 2
rse[, c("A", "D", "F")]
