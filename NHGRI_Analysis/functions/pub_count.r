
install.packages("RISmed")
library(RISmed)

# Counts the number of publications using 2 query words.
pub_count <- function(q1,q2) {
  q <-  EUtilsSummary(paste(q1,q2,sep=" "), type="esearch", db="pubmed")
  print(paste(q1,q2,sep=" "))
  QueryCount(q)
}
