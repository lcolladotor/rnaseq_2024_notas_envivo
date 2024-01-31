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


## Explora el objeto rse de forma interactiva
library("iSEE")
iSEE::iSEE(rse)






initial <- list()

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "gene_101", Search = "", SearchColumns = "",
    HiddenColumns = character(0), VersionInfo = list(iSEE = structure(list(
        c(2L, 14L, 0L)), class = c("package_version", "numeric_version"
    ))), PanelId = c(RowDataTable = 1L), PanelHeight = 500L,
    PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data",
    XAxisColumnData = "Treatment", XAxisFeatureName = "gene_1",
    XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
    YAxisFeatureName = "gene_101", YAxisFeatureSource = "RowDataTable1",
    YAxisFeatureDynamicSource = FALSE, FacetRowByColData = "Treatment",
    FacetColumnByColData = "Treatment", ColorByColumnData = "Treatment",
    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "Treatment", SizeByColumnData = NA_character_,
    TooltipColumnData = character(0), FacetRowBy = "None", FacetColumnBy = "None",
    ColorBy = "Column data", ColorByDefaultColor = "#000000",
    ColorByFeatureName = "gene_1", ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "A",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
    PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    CustomLabels = FALSE, CustomLabelsText = "A", FontSize = 1,
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
    LabelCenters = FALSE, LabelCentersBy = "Treatment", LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(FeatureAssayPlot = 1L),
    PanelHeight = 500L, PanelWidth = 8L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list())

################################################################################
# Settings for Sample assay plot 1
################################################################################

initial[["SampleAssayPlot1"]] <- new("SampleAssayPlot", Assay = "logcounts", XAxis = "None", XAxisRowData = "feature_id",
    XAxisSampleName = "A", XAxisSampleSource = "---", XAxisSampleDynamicSource = FALSE,
    YAxisSampleName = "C", YAxisSampleSource = "ColumnDataTable1",
    YAxisSampleDynamicSource = FALSE, FacetRowByRowData = NA_character_,
    FacetColumnByRowData = NA_character_, ColorByRowData = "feature_id",
    ColorBySampleNameAssay = "logcounts", ColorByFeatureNameColor = "#FF0000",
    ShapeByRowData = NA_character_, SizeByRowData = NA_character_,
    TooltipRowData = character(0), FacetRowBy = "None", FacetColumnBy = "None",
    ColorBy = "None", ColorByDefaultColor = "#000000", ColorByFeatureName = "gene_1",
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "A", ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
    PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    CustomLabels = FALSE, CustomLabelsText = "gene_1", FontSize = 1,
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
    LabelCenters = FALSE, LabelCentersBy = NA_character_, LabelCentersColor = "black",
    VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(SampleAssayPlot = 1L),
    PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list())

################################################################################
# Settings for Column data table 1
################################################################################

initial[["ColumnDataTable1"]] <- new("ColumnDataTable", Selected = "C", Search = "", SearchColumns = "",
    HiddenColumns = character(0), VersionInfo = list(iSEE = structure(list(
        c(2L, 14L, 0L)), class = c("package_version", "numeric_version"
    ))), PanelId = c(ColumnDataTable = 1L), PanelHeight = 500L,
    PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE,
    CustomRowsText = "gene_1\ngene_101\n# gene_1000\ngene_150\ngene_200",
    ClusterRows = TRUE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2",
    DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = "Treatment",
    RowData = character(0), CustomBounds = FALSE, LowerBound = NA_real_,
    UpperBound = NA_real_, AssayCenterRows = TRUE, AssayScaleRows = TRUE,
    DivergentColormap = "purple < black < yellow", ShowDimNames = "Rows",
    LegendPosition = "Bottom", LegendDirection = "Horizontal",
    VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10,
    ShowColumnSelection = TRUE, OrderColumnSelection = TRUE,
    VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version",
    "numeric_version"))), PanelId = c(ComplexHeatmapPlot = 1L),
    PanelHeight = 500L, PanelWidth = 4L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list())

iSEE::iSEE(rse, initial = initial, appTitle = "iSEE pre configurada")


## Descarguemos unos datos de spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")
## Revisemos el tamaño de este objeto
lobstr::obj_size(sce_layer)
iSEE::iSEE(sce_layer)
