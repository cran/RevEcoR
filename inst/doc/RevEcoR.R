## ----setup, include=FALSE------------------------------------------------
library(knitr)
library(RevEcoR)
opts_chunk$set(fig.width=8, fig.height=5)
knit_hooks$set(htmlcap = function(before, options, envir) {
  if(!before) {
    paste('<p class="caption" style="text-align: center; font-size: 20px; color: blue">',options$htmlcap,"</p>",sep="")
    }
    })
set.seed(60823316) 

## ----eval=FALSE----------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite(c("mmnet","KEGGREST","Biobase"))

## ----eval=FALSE----------------------------------------------------------
#  install.packages("RevEcoR")
#  ## or you can install the latest version from github
#  if (!require(devtools)
#    install.packages("devtools")
#  devtools::install_github("yiluheihei/RevEcoR")

## ----eval=TRUE-----------------------------------------------------------
library(RevEcoR) 

## ----eval=FALSE----------------------------------------------------------
#  ## download sample metabolic data from remote KEGG database
#  buc <- getOrgMetabolicData("buc")
#  data(kegg_buc)
#  head(buc)

## ----eval=TRUE, htmlcap="Figure 1 Reconstruction metabolic network of *Buchnera aphidicola APS*", fig.lp="Figure 1", fig.width=8, fig.height=8----
## species in KEGG 
buc.net <- reconstructGsMN(kegg_buc, RefData = NULL) 
igraph::print.igraph(buc.net) 
igraph::plot.igraph(buc.net, vertex.label=NA, vertex.size=5, edge.arrow.size=0.1)
  
## ko annotation profile species detected in a human microbiome in IMG (not in KEGG) 
annodir <- system.file("extdata/koanno.tab",package = "RevEcoR") 
metabolic.data <- read.delim(annodir,stringsAsFactors=FALSE) 
##load the reference metabolic data 
data(RefDbcache,package="mmnet") 
g2 <- reconstructGsMN(metabolic.data, RefData = RefDbcache) 

## ----eval=TRUE, htmlcap="Figure 2The node colored with red represents the species' seed set",fig.lp="Figure 2", fig.width=8, fig.height=8----
## seed set prediction
seed.set <- getSeedSets(buc.net, 0.2) 
show(seed.set) 
head(seed.set@seeds)
## The node colored with red represents the species' seed set
nodes  <- igraph::V(buc.net)$name
seeds  <- unlist(seed.set@seeds)
seed.index  <- match(seeds,nodes)
node.color <- rep("SkyBlue2",length(nodes))
node.color[seed.index]  <- "red"
igraph::plot.igraph(buc.net, 
          vertex.label=NA, vertex.size=5, edge.arrow.size=0.1,
          vertex.color = node.color)

## ------------------------------------------------------------------------
# ptr metabolic network 
data(kegg_ptr) 
##ptr.net <- reconstructGsMN(getOrgMetabolicData("ptr")) 
ptr.net <- reconstructGsMN(kegg_ptr) 
# cooperation analysis between buc and ptr 
cooperation.index <- caculateCooperationIndex(buc.net,ptr.net) 
cooperation.index 

## ------------------------------------------------------------------------
data(gut_microbiome) 
## summary(gut_microbiome) 

## ---- eval = FALSE, echo = TRUE------------------------------------------
#  gut.nets <- lapply(gut_microbiome,reconstructGsMN)
#  seed.sets <- lapply(gut.nets,getSeedSets)
#  gut.interactions <- caculateCooperationIndex(gut.nets)
#  competition.index <- gut.interactions$competition.index
#  complementarity.index <- gut.interactions$complementarity.index

## ---- eval = TRUE, echo = TRUE-------------------------------------------
occurrence.score <- read.delim(system.file("extdata/occurrence.tab",
  package = "RevEcoR"),stringsAsFactors = FALSE, quote = "")

## ---- eval=FALSE,echo=TRUE-----------------------------------------------
#  competition.index <- (competition.index + t(competition.index))/2
#  complementarity.index <- (complementarity.index + t(complementarity.index))/2

## ---- eval=FALSE,echo=TRUE-----------------------------------------------
#  ## upper triangles, which is used to calculate the correlation
#  competition.upper <- competition.index[upper.tri(competition.index)]
#  occurrence.upper <- occurrence.score[upper.tri(occurrence.score)]
#  complementarity.upper <- complementarity.index[upper.tri(complementarity.index)]
#  
#  ## calculate the spearman correlation betwwen co-occurrence scores and two
#  ## interactions indices
#  competition.cor <- cor(competition.upper,occurance.upper,method="spearman")
#  complematarity.cor <- cor(complementarity.upper,occurance.upper,method="spearman")
#  
#  ## permutation-based mantel test. Random permutation the co-occurance score
#  ## 10000 times, P value is the fraction of correlations as high as or higher
#  ## than the original
#  if (require(magrittr)){
#    null.stat <- replicate(10000,
#      sample(1:116) %>% occurrence.score[.,.] %>%
#        .[upper.tri(.)]
#    )
#    competition.null <- cor(competition.upper,null.stat)
#    complementarity.null <- cor(complementarity.upper,null.stat)
#    length(which(competition.null >= competition.cor)) ## 0 p.competition < 0.00001
#    length(which(competition.null >= complementarity.cor)) ## 0 p.complementarity< 0.00001
#  }

## ---- eval=TRUE----------------------------------------------------------
sessionInfo() 

