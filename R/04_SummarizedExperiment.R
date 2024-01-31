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
    assays = SimpleList(
        counts = counts,
        logcounts = log2(counts + 0.5)
    ),
    rowRanges = rowRanges,
    colData = colData
)

## Exploremos el objeto resultante
rse

## Algunos ejemplos de comandos para explorar nuestro
## objeto de RangedSummarizedExperiment (RSE)
dim(rse)
assayNames(rse) # identico a names(assays(rse))
head(assay(rse))
rowRanges(rse)
mcols(rowRanges(rse))
rowData(rse)
rowData(rse)$feature_id[3]
colData(rse)
colData(rse)$Treatment[4]
rse$Treatment[4]

head(assay(rse, "logcounts"))
head(assay(rse, "logcounts")["gene_1", "D"])

## Usamos información de
## https://lcolladotor.github.io/jhustatcomputing2023/posts/14-r-nuts-and-bolts/
l <- list(
    1,
    "a"
)
class(l)
class(l[[1]])
class(l[[2]])

v <- c(
    1,
    "a"
)
class(v)
class(v[1])


l2 <- list(
    "hola" = c(1:3),
    "adios" = letters[1:3]
)
data.frame(
    "hola" = c(1:3),
    "adios" = letters[1:3]
)

l3 <- list(
    "hola" = c(1:4),
    "adios" = letters[1:3]
)
data.frame(
    "hola" = c(1:4),
    "adios" = letters[1:3]
)

class(assays(rse))

assays(rse)[["counts"]]

## Comando 1
rse[1:2, ]
rownames(rse)[1:2]
rse[c("gene_1", "gene_2"), ]
## Comando 2
rse[, c("A", "D", "F")]
rse[, c(1, 4, 6)]
colnames(rse)[c(1, 4, 6)]

## Uso de memoria
rse2 <- rse[, c(1, 4, 6)]
ls()
rm(rse2)
?gc
