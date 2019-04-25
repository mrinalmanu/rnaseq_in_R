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
#######

# SOLUTION: the new script should look something like this

install.packages("plyr")
library(plyr)


ddply(.data, .variables, .fun = NULL, ..., .progress =
  "none",.inform = FALSE, .drop = TRUE, .parallel = FALSE,
  .paropts = NULL)

# They key is the first three arguments. 
# The first one is the dataframe you want to  work with. 
# The second one is the column name to get rows with the same information. 
# In the example above, that would be “XXXXXX”. 
# The third argument, function, is the most critical one. 
# The command will first look at the top cell of the column specified by the second parameter.
# Then it will find all rows with the same information. 
# For example, assume the following  ddply function is called

avgdata<-ddply(data1,.(entrezID),calavg)

# If data1 is the dataframe of the example on the top of page,
# ddply will grab the first entry in the second argument “entreZID”,
# which is “24244”. Then it will look for other rows that also contain
# “24244” in entrezID. In this example, there are two rows with the same “entrezID”, 
# it will send a dataframe with theses two rows to the function of the third argument. 

calavg<-function(x)
   avg<-apply(x[,2:3],2,mean)
   return(avg)
}
# This simple function takes a dataframe and applies column-wide  mean function.
# The results are stored in a dataframe “avg”, which will be returned upon completion. 

# https://bioinfomagician.wordpress.com/2014/07/27/my-favorite-commands-part2-data-reduction-using-ddply/

