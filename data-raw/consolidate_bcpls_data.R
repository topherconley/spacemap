
## Distill BCPLS Analysis into a small list of objects 
## ready to be input to network analysis toolkit of spacemap. 



####(1)
#prepare annotation 

datadir <- "~/repos/neta-bcpls/data/"
library(Biobase)
##X information
cnaset <- readRDS(file = file.path(datadir,"cna-expression-set.rds"))
xy_node_attributes <- featureData(cnaset)
xy_node_attributes <- pData(xy_node_attributes)
xinfo <- xy_node_attributes[,c(1,9,2,3,4)]
xinfo$start <- as.integer(xinfo$start)
xinfo$end <- as.integer(xinfo$end)
names(xinfo)[1] <- "id"
names(xinfo)[2] <- "alias"
rownames(xinfo) <- NULL
head(xinfo)


##Y information
protset <- readRDS(file = file.path(datadir, "prot-expression-set.rds"))
yy_node_attributes <- featureData(protset)
yy_node_attributes <- pData(yy_node_attributes)
yy_node_attributes$description <- as.character(yy_node_attributes$description)
yinfo <- yy_node_attributes[,c(1,7,2,3,4,8,10)]
yinfo$start <- as.integer(yinfo$start)
yinfo$end <- as.integer(yinfo$end)
names(yinfo)[1] <- "id"
names(yinfo)[2] <- "alias"
rownames(yinfo) <- NULL
head(yinfo)

####(2)
#network 
moddir <- "~/repos/neta-bcpls/model-fits/"
suppressMessages(library(Matrix))
net <- readRDS(file = file.path(moddir, "smap_prot_boot_vote.rds"))
names(net) <- c("yy", "xy")
rownames(net$xy) <- xinfo$id
colnames(net$xy) <- yinfo$id
rownames(net$yy) <- yinfo$id
colnames(net$yy) <- yinfo$id
str(net)
#net$yy <- as.matrix(net$yy)
#net$xy <- as.matrix(net$xy)

###(3)
#bootstrap degree distributions
bdeg <- readRDS(file = file.path(moddir, "smap_prot_boot_degree.rds"))
str(bdeg)
names(bdeg) <- c("xy", "yy")
dim(bdeg$yy)
colnames(bdeg$xy) <- xinfo$id
colnames(bdeg$yy) <- yinfo$id

str(bdeg)
###(4)
#GO 2 entrez ID
go2eg <- readRDS(file = file.path(datadir, 
                                  "prot-go-bp-to-entrez-gene-list.rds"))

bcpls <- list(net = net, 
              yinfo = yinfo,
              xinfo = xinfo,
              bdeg =  bdeg, 
              go2eg = go2eg)

library(devtools)
devtools::use_data(bcpls, overwrite = T)

#check that the total .rda data objects is less than 1MB.
stopifnot(sum(tools::checkRdaFiles("data/")$size) < 1e6)
