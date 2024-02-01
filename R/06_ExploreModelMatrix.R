## ?model.matrix
mat <-
    with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat

class(mat)
dim(mat)
dim(trees)
head(trees)
head(mat)


head(model.matrix(log(trees$Volume) ~ log(trees$Height) + log(trees$Girth)))

summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))

## Ejemplo 1 de ExploreModelMatrix
library("ExploreModelMatrix")
(sampleData <- data.frame(genotype = rep(c("A", "B"), each = 4),
    treatment = rep(c("ctrl", "trt"), 4)))

app <- ExploreModelMatrix(sampleData = sampleData,
    designFormula = ~ genotype + treatment)
if (interactive())
    shiny::runApp(app)


## Ejemplo 2 de ExploreModelMatrix
(
    sampleData <- data.frame(
        Response = rep(c("Resistant", "Sensitive"), c(12, 18)),
        Patient = factor(rep(c(1:6, 8, 11:18), each = 2)),
        Treatment = factor(rep(c("pre", "post"), 15)),
        ind.n = factor(rep(c(1:6, 2, 5:12), each = 2))
    )
)

app <- ExploreModelMatrix(
    sampleData = sampleData,
    designFormula = ~ Response + Response:ind.n + Response:Treatment
)
if (interactive())
    shiny::runApp(app)

## Ejemplo 3 de ExploreModelMatrix
(sampleData = data.frame(condition = factor(rep(
    c("ctrl_minus", "ctrl_plus",
        "ko_minus", "ko_plus"), 3
)),
    batch = factor(rep(1:6, each = 2))))

app <- ExploreModelMatrix(sampleData = sampleData,
    designFormula = ~ batch + condition)

if (interactive())
    shiny::runApp(app)
