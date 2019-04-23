
# read the dataset into R
library(GEOquery)
library(limma)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

# for collapseBy
source("//home/manu/Documents/assignment/dataset_3_neuro/functions.r")

dir.create("cache")

res <- getGEO("GSE28790", AnnotGPL = TRUE, destdir="cache")[[1]]
# info: SIRT1 impact on global gene expression in the brain


###Collapsing the data
# for collapseBy
source("//home/student/deseq/functions.R")

str(experimentData(res))
str(pData(res))
head(fData(res))
res$`title`

res$condition <- c("BSKO", "BSKO", "BSKO", "BSKO",
                   "WT", "WT", "WT", "WT")
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

fit <- lmFit(res.qnorm.top12K, res.design)

fit2 <- contrasts.fit(fit, makeContrasts(conditionBSKO-conditionWT, 
                                         levels=res.design))

fit2 <- eBayes(fit2)
de <- topTable(fit2, adjust.method="BH", number=Inf)
head(de)

library(data.table)
de <- as.data.table(de, keep.rownames=TRUE)
de[entrez == "Sirt2"]
de[entrez == "Sirt3"]
de[entrez == "Sirt4"]
de[entrez == "Sirt7"]


# ####
# FGSEA
# ####

de2 <- data.frame(de$entrez, de$t)
colnames(de2) <- c('ENTREZ', 'stat')

library(fgsea)

ranks <- deframe(de2)
head(ranks, 20)

# Load the pathways into a named list
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

pdf("fgseaResTidy.pdf", width = 20, height = 20)

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

dev.off()


# We can see that two of the hallmark pathways are significant

# let's look at all pathways and select for 
# pathways associated with SIRT in significance threshold
# running the fgsea algorithm on all pathways

fgseaRes.all <- fgsea(pathways=pathways, stats=ranks, nperm=1000)

# number <- data.frame(grep("Sirt", fgseaRes.all$leadingEdge))
# colnames(number) <- c('row_number')
# Sirt <- subset(fgseaRes.all, rownames(fgseaRes.all) %in% number$row_number)

# using tidy to view pretty results

fgseaResTidy <- fgseaRes.all %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table for all pathways
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

# ggplotting for all pathways

pdf("fgseaResTidy_all.pdf", width=20, height=500)

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="All pathways NES from GSEA") + 
  theme_minimal()

dev.off()

# plot pca

df_12k <- data.frame(res.qnorm.top12K@assayData$exprs)
colnames(df_12k) <- c('BSKO', 'BSKO', 'BSKO', 'BSKO',
                      'WT', 'WT', 'WT', 'WT')
p <- prcomp(na.omit(df_12k))

devtools::install_github("sinhrks/ggfortify")
library(ggfortify)

ggplot2::autoplot(p, label = FALSE, shape = FALSE, 
                  loadings.label = TRUE)


ggplot(stack(df_12k), aes(x = ind, y = values)) +
  geom_boxplot()

# We can see that there are various activated pathways

# # SIRT1 is a NAD dependent deacetylase that provides coping 
# mechanism against nutritional changes in cell. 
# It is shown that by activation of genes encoding for 
# MAO-A (monoamine oxidase). SSRIs are inhibited. 
# SSRIs play a key role in uptake of serotonin and 
# thus manifestation of symptoms of depression.
#  
# # In analysis we can see that indeed expression 
# of this gene is increased. Further, MAO-A inhibitors or SSRIs 
# (selective serotonin reuptake inhibitors) 
# normalise anxiety differences between wild-type and mutant animals.
# # From network analysis we can see the involvement 
# in Sleep Cycle. Perhaps a disrupted sleep cycle and 
# altered cell signalling response is promoting 
# anxious behaviour in these mice.
# 
