library(RMySQL)
# Define the function needed for this script
# Find table of analyte has pathway from given pathway IDs
# Aggregate ramp Id to ramp pathway Id 
# GC is C or G
findAnalyteHasPathway <- function(pathwayRampId,GC = "C",n = 10){
  con <- dbConnect(MySQL(),
                   user = 'root',
                   dbname='mathelabramp',
                   password = 'Ehe131224',
                   host = 'localhost')
  on.exit(dbDisconnect(con))
  p_id <- unique(pathwayRampId)
  p_id <- sapply(p_id,shQuote)
  p_id <- paste(p_id,collapse = ",")
  query <-paste0("select * from analytehaspathway where pathwayRampId in (",
                 p_id,
                 ");")
  df <- dbGetQuery(con,
                   query)
  df2 <- aggregate(df$rampId,list(df$pathwayRampId),FUN = function(x){
    x <- x[grepl(paste0("RAMP_",GC,"_"),x)]
    if(length(x) >= n ){
      paste(x,collapse = ",")
    } else {
      x <- 0
    }
  })
  fdf <- df2[df2$x!=0,]
  fdf2 <- data.frame(fdf[,-1],row.names = fdf[,1],stringsAsFactors = F)
  df.list <- setNames(split(fdf2, seq(nrow(fdf2))), rownames(fdf2))
  df.list <- lapply(df.list,FUN = function(x){
    text <- x[[1]]
    text <- strsplit(text,split = ",")
  })
  df.list <- lapply(df.list,unlist)
} # end of function

# Build connection with the database
# user, dbname, password, host depdent on what you have
con <- dbConnect(MySQL(),
                 user = 'root',
                 dbname='mathelabramp',
                 password = 'Ehe131224',
                 host = 'localhost')

pathways<- dbGetQuery(con,'select * from pathway;')


dbname <- unique(pathways$type)

# pathwayInHmdb <- pathways[pathways$type == 'hmdb',]
pathwayInKegg <- pathways[pathways$type == 'kegg',]
pathwayInWiki <- pathways[pathways$type == 'wiki',]
pathwayInReac <- pathways[pathways$type == 'reactome',]

# define the minimum metabolites/genes that the pathway shuold have to take it in account
# The paper use 5.
# The number is higher, the below part takes shorter time.
min_analyte <- 5

# Store Compound Ids in List
# listOfHmdbC <- findAnalyteHasPathway(pathwayInHmdb$pathwayRampId)
listOfKeggC <- findAnalyteHasPathway(pathwayInKegg$pathwayRampId,n = min_analyte)
listOfWikiC <- findAnalyteHasPathway(pathwayInWiki$pathwayRampId,n = min_analyte)
listOfReacC <- findAnalyteHasPathway(pathwayInReac$pathwayRampId,n = min_analyte)
# Store Gene Ids in List
# listOfHmdbG <- findAnalyteHasPathway(pathwayInHmdb$pathwayRampId,GC="G")
listOfKeggG <- findAnalyteHasPathway(pathwayInKegg$pathwayRampId,GC="G",n = min_analyte)
listOfWikiG <- findAnalyteHasPathway(pathwayInWiki$pathwayRampId,GC="G",n = min_analyte)
listOfReacG <- findAnalyteHasPathway(pathwayInReac$pathwayRampId,GC="G",n = min_analyte)
# Setup minimum number of analytes that will be considered

# May need to filter out that pathway that has less than 5 metabolites 
# Output to a matrix
# In order of HMDB Kegg Wiki Reac

pathwayid <- c(#names(listOfHmdbC),
  names(listOfKeggC),
  names(listOfWikiC),
  names(listOfReacC))


metabolite_result <- matrix(NA,nrow = length(pathwayid),ncol = length(pathwayid))

# Assign names on the metabolites result
colnames(metabolite_result) <- pathwayid
rownames(metabolite_result) <- pathwayid
colnames(metabolite_result2) <- pathwayid
rownames(metabolite_result2) <- pathwayid
# Since HMDB SMPDB's behavior is not clear
# We skip the HMDB for the heatmap plot
pathToanalC <- do.call(c,list(#listOfHmdbC,
  listOfKeggC,
  listOfWikiC,
  listOfReacC))
# Compute for the matrix
for(i in 1:length(pathwayid)){
  id <- pathwayid[i]
  cid <- pathToanalC[[i]]
  for (j in 1:length(pathwayid)) {
    if(is.na(metabolite_result[i,j])){
      if(i==j){
        metabolite_result[i,j] <- 1
      }else{
        cid2 <- pathToanalC[[j]]
        shared_metabolite <- unique(intersect(cid,cid2))
        total <- unique(union(cid,cid2))
        metabolite_result[i,j] <- length(shared_metabolite)/length(total)
        print(metabolite_result[i,j])
        if(is.na(metabolite_result[j,i])){
          metabolite_result[j,i] <- metabolite_result[i,j]
        }
      }
    }
    print(paste("Compute for ",i,",",j))
  }
}

### Part for genes ...
## The process is analog to the metabolites part
pathwayidG <- c(#names(listOfHmdbG),
  names(listOfKeggG),
  names(listOfWikiG),
  names(listOfReacG))
gene_result <- matrix(NA,nrow = length(pathwayidG),ncol = length(pathwayidG))

# In order of HMDB Kegg Wiki Reac
# Assign names on the metabolites result
colnames(gene_result) <- pathwayidG
rownames(gene_result) <- pathwayidG
colnames(gene_result2) <- pathwayidG
rownames(gene_result2) <- pathwayidG

pathToanalG <- do.call(c,list(#listOfHmdbG,
  listOfKeggG,
  listOfWikiG,
  listOfReacG))
# Compute for the matrix
for(i in 1:length(pathwayidG)){
  id <- pathwayidG[i]
  cid <- pathToanalG[[i]]
  
  for (j in 1:length(pathwayidG)) {
    if(is.na(gene_result[i,j])){
      if(i==j){
        gene_result[i,j] <- 1
      }else{
        cid2 <- pathToanalG[[j]]
        shared_metabolite <- unique(intersect(cid,cid2))
        total <- unique(union(cid,cid2))
        gene_result[i,j] <- length(shared_metabolite)/length(total)
        print(gene_result[i,j])
        if(is.na(gene_result[j,i])){
          gene_result[j,i] <- gene_result[i,j]
        }
        
      }
    }
    print(paste("Compute for ",i,",",j))
  }
}

library(plotly)
library(magrittr)
# Generate metabolite plot
p_metabolites <- plot_ly(z = metabolite_result,
                         x = pathwayid,
                         y = pathwayid,
                         colors = c('#ece7f2','#2b8cbe'),
                         width = 1200,
                         height = 1200,
                         type = "heatmap") %>%
  layout(title = "metabolites overlap",
         font = list(
           size = 40
         ),
         autosize = TRUE,
         margin = list(l = 150,
                       r = 100,
                       t = 100,
                       b = 200),
         xaxis = list(
           type = "category",
           autoticks = FALSE,
           ticks = "inside",
           showline = FALSE,
           zeroline = FALSE,
           showticklabels = FALSE,
           autorange = TRUE
         ),
         yaxis = list(
           type = "category",
           ticks = "inside",
           showline = FALSE,
           zeroline = FALSE,
           showticklabels =FALSE,
           autorange = TRUE
         )
  )
# Generate gene plot.
p_genes <- plot_ly(z = gene_result,
                   x = pathwayidG,
                   y = pathwayidG,
                   width = 2000,
                   height = 2000,
                   colors = c('#ece7f2','#2b8cbe'),
                   type = "heatmap") %>%
  layout(title = "Genes overlap",
         font = list(
           size = 40
         ),
         autosize = TRUE,
         margin = list(l = 150,
                       r = 100,
                       t = 100,
                       b = 200),
         xaxis = list(
           type = "category",
           autoticks = FALSE,
           ticks = "inside",
           showline = FALSE,
           zeroline = FALSE,
           showticklabels = FALSE,
           autorange = TRUE
         ),
         yaxis = list(
           type = "category",
           ticks = "inside",
           showline = FALSE,
           zeroline = FALSE,
           showticklabels =FALSE,
           autorange = TRUE
         )
  )