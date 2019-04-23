
# read the dataset into R
library(GEOquery)
library(limma)
#library(org.Mm.eg.db)
library(org.Hs.eg.db)

# for collapseBy
source("//home/manu/Documents/assignment/dataset_3_neuro/functions.r")

dir.create("cache")

res <- getGEO("GSE53890", AnnotGPL = TRUE, destdir="cache")[[1]]
# info: Age effect on normal adult brain: frontal cortical region


###Collapsing the data
# for collapseBy
source("//home/student/deseq/functions.R")

str(experimentData(res))
str(pData(res))
head(fData(res))
res$`Sex:ch1`

res$condition <- gsub("\\+", "_", res$`Sex:ch1`)
res$condition

res <- collapseBy(res, fData(res)$`Gene symbol`, FUN=median)
res <- res[!grepl("///", rownames(res)), ]
res <- res[rownames(res) != "", ]

# there is a lot of garbage there



fData(res) <- data.frame(row.names = rownames(res))

fData(res)$entrez <- row.names(fData(res))

fData(res)$symbol <- mapIds(org.Hs.eg.db, keys=fData(res)$entrez, 
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

fit2 <- contrasts.fit(fit, makeContrasts(conditionFemale-conditionMale, 
                                         levels=res.design))

fit2 <- eBayes(fit2)
de <- topTable(fit2, adjust.method="BH", number=Inf)
head(de)
library(data.table)
de <- as.data.table(de, keep.rownames=TRUE)
de[entrez == "REST"]




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

pdf("fgseaResTidy.pdf", width = 100, height = 1000)

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

dev.off()



# let's look at all pathways and select for 
# pathways associated with brain in significance threshold
# running the fgsea algorithm on all pathways

fgseaRes.all <- fgsea(pathways=pathways, stats=ranks, nperm=1000)

number <- data.frame(grep("REST", fgseaRes.all$leadingEdge))
colnames(number) <- c('row_number')
REST <- subset(fgseaRes.all, rownames(fgseaRes.all) %in% number$row_number)

# using tidy to view pretty results

fgseaResTidy <- REST %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table for all pathways
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

# ggplotting for all pathways

pdf("fgseaResTidy_REST.pdf", width=100, height=1000)

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="All pathways NES from GSEA") + 
  theme_minimal()

dev.off()

# We can see that two pathways are significant, associated with 
# ion channel activity

#scale rows
xt<-t(as.matrix(res.qnorm.top12K))
xts<-scale(xt)
xtst<-t(xts)
xtst <- na.omit(xtst)
colnames(xtst) <- res$`Sex:ch1`

#only grab top 1000 by p-value
h<-head(xtst, n = 1000L)

#set layout options - adjust if labels get cut off
pdf("heatmap.pdf",width=500, height=500)

#draw heatmap allowing larger margins and adjusting row label font size
heatmap(h)

#output plot to file
dev.off()

# install.packages('devtools')
library(devtools)
# devtools::install_github("sinhrks/ggfortify")
library(ggfortify)

pdf('box_dataset1.pdf')

ggplot(stack(data.frame(t(xt))), aes(x = ind, y = values)) +
  geom_boxplot()

dev.off()
