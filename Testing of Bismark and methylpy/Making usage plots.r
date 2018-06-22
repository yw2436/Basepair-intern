##
get_first <- function(tmp){
  tmp <- strsplit(tmp, split = ' ')
  tmp <- tmp[[1]][tmp[[1]]!=''][1]
  as.integer(tmp)
}
##

me_cpu_10 <- read.delim('methylpy_cpu_usage_10X.txt', sep = ':', header = FALSE,as.is = TRUE)
me_cpu_10 <- me_cpu_10[-c(1,2,3),]
me_cpu_10 <- apply(me_cpu_10,MARGIN = c(1,2), get_first)


me_cpu_30 <- read.delim('methylpy_cpu_usage_30X.txt', sep = ':', header = FALSE,as.is = TRUE)
me_cpu_30 <- me_cpu_30[-c(1,2,3),]
me_cpu_30 <- apply(me_cpu_30,MARGIN = c(1,2), get_first)

##
get_four <- function(x){
  x[x!='' & !is.na(x)]
}
##

me_mem_10 <- read.delim('methylpy_mem_usage_10X.txt', skip = 1, header = TRUE, sep = ' ', as.is = T)
me_mem_10 <- t(apply(me_mem_10, 1, get_four))
colnames(me_mem_10) <- c('used','buff','cach','free')

me_mem_30 <- read.delim('methylpy_mem_usage_30X.txt', skip = 1, header = TRUE, sep = ' ', as.is = T)
me_mem_30 <- t(apply(me_mem_30, 1, get_four))
colnames(me_mem_30) <- c('used','buff','cach','free')

