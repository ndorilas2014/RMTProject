library(tidyverse)

load("./dataestimate_time.Rdata")

p <- proj_time %>% select(sname, experiment_day, otu_id, count)

for ( id in unique(p$sname) ){
  idd <- gsub(" ", "_", id)
  dmat <- p %>% filter(sname == id) %>% spread( otu_id, count, fill = 0 ) %>% select(-sname) %>% 
    column_to_rownames('experiment_day') %>% as.matrix() %>% t()
  assign(idd, dmat)  
  save( file =  paste("./",idd,".Rdata", sep = ""), list=idd )  
}