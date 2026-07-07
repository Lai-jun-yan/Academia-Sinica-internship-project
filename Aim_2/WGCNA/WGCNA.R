set.seed(123)
# sample size
nSamples <- 100

# module eigengene
eigengenes <- data.frame(
  MEblue = rnorm(nSamples),
  MEbrown = rnorm(nSamples),
  MEgreen = rnorm(nSamples),
  MEyellow = rnorm(nSamples),
  MEred = rnorm(nSamples)
)

modProportions <- c(
  0.2,
  0.2,
  0.2,
  0.2,
  0.1,
  0.1
)

sim <- simulateDatExpr(
  eigengenes = eigengenes,
  nGenes = 1000,
  modProportions = modProportions
)

expr <- sim$datExpr

powers <- 1:20

sft <- pickSoftThreshold(
  expr,
  powerVector=powers
)

net <- blockwiseModules(
  expr,
  power=6,
  minModuleSize=30
)

barplot(
  table(sim$setLabels),
  main="True simulated modules",
  xlab="Module",
  ylab="Number of genes"
)

corMat <- cor(expr)
corMat[1:5,1:5]

beta <- 6

adjacency <- adjacency(
  expr,
  power=beta
)

TOM <- TOMsimilarity(
  adjacency
)
TOM[1:5,1:5]

dissTOM <- 1-TOM

geneTree <- hclust(
  as.dist(dissTOM),
  method="average"
)

plot(
  geneTree,
  main="Gene clustering"
)

library(dynamicTreeCut)

modules <- cutreeDynamic(
  dendro=geneTree,
  distM=dissTOM,
  deepSplit=2,
  minClusterSize=30
)

table(modules)

colors <- labels2colors(modules)
MEs <- moduleEigengenes(
  expr,
  colors
)$eigengenes

head(MEs)

trait <- data.frame(
  Disease=rbinom(
    100,
    1,
    0.5
  )
)

moduleTraitCor <- cor(
  MEs,
  trait
)

labeledHeatmap(
  Matrix=moduleTraitCor,
  xLabels=names(trait),
  yLabels=names(MEs)
)





