
# read the dataset into R
library(GEOquery)
library(limma)
# library(org.Mm.eg.db)
library(org.Hs.eg.db)

# for collapseBy
source("//home/manu/Documents/assignment/dataset_3_neuro/functions.r")

dir.create("cache")

# to get the dataset uncomment the comment line below and execute

res <- getGEO("GSE19404", AnnotGPL = TRUE, destdir="cache")[[1]]
# info: SIRT1 impact on global gene expression in the brain
# 

###Collapsing the data
# for collapseBy

str(experimentData(res))
str(pData(res))
head(fData(res))
res$characteristics_ch1 <- NULL
res$characteristics_ch1 <- t(annotate[1,])

res$`title`

annotate <- t(annotate)
res$condition <- NULL
res$condition <- annotate
res$condition


res <- collapseBy(res, fData(res)$`Gene symbol`, FUN=median)
res <- res[!grepl("///", rownames(res)), ]
res <- res[rownames(res) != "", ]


# there is a lot of garbage there

fData(res) <- data.frame(row.names = rownames(res))

fData(res)$entrez <- row.names(fData(res))

fData(res)$symbol <- mapIds(org.Mm.eg.db, keys=fData(res)$entrez, 
                            keytype="SYMBOL", column="ENTREZID" )

res.qnorm <- res

summary(exprs(res.qnorm))
exprs(res.qnorm) <- normalizeBetweenArrays(log2(exprs(res.qnorm)+1), method="quantile")
summary(exprs(res.qnorm))

res.qnorm.top12K <- res.qnorm
res.qnorm.top12K <- res.qnorm.top12K[head(order(apply(exprs(res.qnorm.top12K), 1, mean), 
                                                decreasing = TRUE), 12000), ]

res.design <- model.matrix(~0+condition, data=pData(res.qnorm.top12K))

intermediate <- data.frame(res.design)
colnames(intermediate) <- c("ATRT", "CNS", "Medulloblastoma", "Normal", "Pineoblastoma")
rm(res.design)
res.design <- as.matrix(intermediate)

fit <- lmFit(res.qnorm.top12K, res.design)

#############
fit2 <- contrasts.fit(fit, makeContrasts(ATRT-Normal,
                                         CNS-Normal,
                                         Medulloblastoma-Normal,
                                         Pineoblastoma-Normal,
                                         levels=res.design))

fit2 <- eBayes(fit2)
de <- topTable(fit2, adjust.method="BH", number=Inf)
head(de)

library(data.table)
de <- as.data.table(de, keep.rownames=TRUE)
entry <- function(){
  item <- data.frame("CDK4","SEC61G","TSPAN31","LANCL2",
                     "EGFR","APC","ATM","BMPR1A","BRCA1",
                     "BRCA2","CDK4","CDKN2A","CREBBP","EGFR",
                     "EP300","ETV6","FHIT","FLT3","HRAS","KIT",
                     "MET","MLH1","NTRK1","PAX8","PDGFRA",
                     "PRCC","PRKAR1A","PTEN","RET","STK11",
                     "TFE3","TP53","WWOX")
  item<- t(item)
  rownames(item) <- NULL
  
  x<- for (i in item){ 
    print(de[entrez == i])
    
  }
  
  return(x)
  
}

entry()

# ####
# FGSEA
# ####
# install.packages('fgsea')
# install.packages('tibble')

library(fgsea)
library(tibble)

de2 <- data.frame(de$entrez, de$P.Value)
colnames(de2) <- c('ENTREZ', 'stat')

ranks <- deframe(de2)
head(ranks, 20)

# Load the pathways into a named list

# install.packages('msigdbr')
library(msigdbr)

m_df <- msigdbr(species = "Homo sapiens")
m_df
pathways <- split(m_df$human_gene_symbol, m_df$gs_name)

# filter the list to include only hallmark pathways

library(data.table)

pathways.hallmark <- m_df[m_df$gs_name %like% "HALLMARK_", ]
pathways.hallmark <- split(pathways.hallmark$human_gene_symbol, pathways.hallmark$gs_name)

# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark %>% 
  head() %>% 
  lapply(head)


# running the fgsea algorithm on hallmark.pathways

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
# ggplotting for hallmark pathways

library(ggplot2)
pdf("fgseaResTidy.pdf", width = 20, height = 20)

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

dev.off()


# We can sees seven significant hallmark pathways

# let's look at all pathways and select for 
# pathways associated with genes of intrest in significance threshold
# running the fgsea algorithm on all pathways

fgseaRes.all <- fgsea(pathways=pathways, stats=ranks, nperm=1000)

# searching for the genes in pathway and appending the rownumbers
sink('numbers.txt')

options(max.print=2000)
for(i in item){
  print(grep(i, fgseaRes.all$leadingEdge))
}

sink()

# we have to do a lot of cleaning of the data before importing it as csv
# getting only unique values from all numbers
unique_vals <- data.frame(as.integer(unique(unlist(numbers))))

colnames(unique_vals) <- c('row_number')
final <- subset(fgseaRes.all, rownames(fgseaRes.all) %in% unique_vals$row_number)

# using tidy to view pretty results
# install.packages('DT')

fgseaResTidy <- fgseaRes.all %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table for all pathways
library(DT)

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()


# ggplotting for all pathways

pdf("fgseaResTidy_all.pdf", width=20, height=2000)

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="All pathways NES from GSEA") + 
  theme_minimal()

dev.off()

# plot pca

df_12k <- data.frame(res.qnorm.top12K@assayData$exprs)
colnames(df_12k) <- annotate
p <- prcomp(na.omit(df_12k))
 
# install.packages('devtools')
library(devtools)
# devtools::install_github("sinhrks/ggfortify")
library(ggfortify)

pdf('pca_and_box_dataset3.pdf')

ggplot2::autoplot(p, label = TRUE, shape = FALSE, 
                  loadings.label = TRUE)

ggplot(stack(df_12k), aes(x = ind, y = values)) +
  geom_boxplot()

dev.off()
#
# In the paper, authors identified similarities 
# in transcriptomes of distinct human brain tumor, 
# specifically glioblastomas of TCGA classical 
# subtype and oligodendroglial tumors, Philips 
# proferative gliomas and AT/RT.