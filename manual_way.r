# Manual script

library(dplyr)

new_GSE <- data.frame(GSE19404_series_matrix)
new_GSE$ID_REF = NULL

res <- new_GSE
res$id <- meta_data$`Gene ID`
# res$symbol <- meta_data$`Gene symbol`

es <- res 

# removing the annoying '///' character that is problematic in es$id column
es$id <- gsub('///', '', es$id)

# converting the class of es$id from character to int
es$id <- as.integer(es$id)

# also dropping the 'NA' values from es$id column
es <- na.omit(es)

# there is some problem after this section
###################################################################
#removing the duplicates and replacing by median of the column
new_es <- es %>% group_by(es$id) %>% summarise_all(funs(median))

# res <- res[rownames(res) != "", ]

new_es$`id` <- NULL

colnames(new_es) <- annotate

#
# Annotation
#
newer_es <- new_es

