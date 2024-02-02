## Código copiado de
## https://github.com/lcolladotor/rnaseq_LCG-UNAM_2024/blob/main/05_modelos.R


## ----download_SRP045638----------------
library("recount3")

human_projects <- available_projects()

rse_gene_SRP045638 <- create_rse(
    subset(
        human_projects,
        project == "SRP045638" & project_type == "data_sources"
    )
)
assay(rse_gene_SRP045638, "counts") <- compute_read_counts(rse_gene_SRP045638)

## Explorar
rse_gene_SRP045638

## Esto no corre
expand_sra_attributes(rse_gene_SRP045638)

rse_gene_SRP045638$sra.sample_attributes[1:3]

## Vamos a eliminar la info de dev_stage
rse_gene_SRP045638$sra.sample_attributes <- gsub("dev_stage;;Fetal\\|", "", rse_gene_SRP045638$sra.sample_attributes)
rse_gene_SRP045638$sra.sample_attributes[1:3]
expand_sra_attributes(rse_gene_SRP045638)

## Si corre, asi que ahora si podemos guardar el
## resultado
rse_gene_SRP045638 <- expand_sra_attributes(rse_gene_SRP045638)

colData(rse_gene_SRP045638)[
    ,
    grepl("^sra_attribute", colnames(colData(rse_gene_SRP045638)))
]

## Pasar de character a numeric o factor
rse_gene_SRP045638$sra_attribute.age <- as.numeric(rse_gene_SRP045638$sra_attribute.age)
rse_gene_SRP045638$sra_attribute.disease <- factor(tolower(rse_gene_SRP045638$sra_attribute.disease))
rse_gene_SRP045638$sra_attribute.RIN <- as.numeric(rse_gene_SRP045638$sra_attribute.RIN)
rse_gene_SRP045638$sra_attribute.sex <- factor(rse_gene_SRP045638$sra_attribute.sex)

## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP045638)[
    ,
    grepl("^sra_attribute.[age|disease|RIN|sex]", colnames(colData(rse_gene_SRP045638)))
]))

## Encontraremos diferencias entre muestra prenatalas vs postnatales
rse_gene_SRP045638$prenatal <-
    factor(ifelse(
        rse_gene_SRP045638$sra_attribute.age < 0,
        "prenatal",
        "postnatal"
    ))
table(rse_gene_SRP045638$prenatal)

## http://rna.recount.bio/docs/quality-check-fields.html
rse_gene_SRP045638$assigned_gene_prop <- rse_gene_SRP045638$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP045638$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP045638$assigned_gene_prop)

with(colData(rse_gene_SRP045638), plot(assigned_gene_prop, sra_attribute.RIN))

## Hm... veamos si hay una diferencia entre los grupos
with(colData(rse_gene_SRP045638), tapply(assigned_gene_prop, prenatal, summary))

## Guardemos nuestro objeto entero por si luego cambiamos de opinión
rse_gene_SRP045638_unfiltered <- rse_gene_SRP045638

## Eliminemos a muestras malas
hist(rse_gene_SRP045638$assigned_gene_prop)

rse_gene_SRP045638 <- rse_gene_SRP045638[, rse_gene_SRP045638$assigned_gene_prop > 0.3]

gene_means <- rowMeans(assay(rse_gene_SRP045638, "counts"))
summary(gene_means)

## Eliminamos genes
rse_gene_SRP045638 <- rse_gene_SRP045638[gene_means > 0.1, ]

## Dimensiones finales
dim(rse_gene_SRP045638)

## Porcentaje de genes que retuvimos
round(nrow(rse_gene_SRP045638) / nrow(rse_gene_SRP045638_unfiltered) * 100, 2)

library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
    counts = assay(rse_gene_SRP045638, "counts"),
    genes = rowData(rse_gene_SRP045638)
)
dge <- calcNormFactors(dge)


mod <- model.matrix(~ prenatal + sra_attribute.RIN + sra_attribute.sex + assigned_gene_prop,
    data = colData(rse_gene_SRP045638)
)
colnames(mod)
head(mod)

library("limma")
vGene <- voom(dge, mod, plot = TRUE) # 2015
eb_results <- eBayes(lmFit(vGene)) # microarreglos

de_results <- topTable(
    eb_results,
    coef = "prenatalprenatal", # 2,
    number = nrow(rse_gene_SRP045638),
    sort.by = "none"
)
dim(de_results)
head(de_results)


## Extraer valores de los genes de interés
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

rownames(exprs_heatmap) <- rowRanges(rse_gene_SRP045638)$gene_name[
    match(
        rownames(exprs_heatmap),
        rowRanges(rse_gene_SRP045638)$gene_id # rownames(rse_gene_SRP045638)
    )
]

## Creemos una tabla con información de las muestras
## y con nombres de columnas más amigables
df <- as.data.frame(colData(rse_gene_SRP045638)[, c("prenatal", "sra_attribute.RIN", "sra_attribute.sex")])
colnames(df) <- c("AgeGroup", "RIN", "Sex")

## Hagamos un heatmap
library("pheatmap")
pheatmap(
    exprs_heatmap,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = df
)

## Versión con centering y scaling en los renglones (los genes)
pheatmap::pheatmap(
    exprs_heatmap,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = df,
    scale = "row"
)


pheatmap::pheatmap(
    t(exprs_heatmap),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    annotation_row = df,
    scale = "column"
)
