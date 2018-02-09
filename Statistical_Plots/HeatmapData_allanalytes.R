library(RMySQL)
# define function ...
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
  if(GC == 'both'){
    df2 <- aggregate(df$rampId,list(df$pathwayRampId),FUN = function(x){
      if(length(x) >= n){
        paste(x,collapse = ',')
      } else{
        x <- 0
      }
    })
  }
  else if (GC %in% c('G','C')){
    df2 <- aggregate(df$rampId,list(df$pathwayRampId),FUN = function(x){
      x <- x[grepl(paste0("RAMP_",GC,"_"),x)]
      if(length(x) >= n ){
        paste(x,collapse = ",")
      } else {
        x <- 0
      }
    })
  }
  fdf <- df2[df2$x!=0,]
  fdf2 <- data.frame(fdf[,-1],row.names = fdf[,1],stringsAsFactors = F)
  df.list <- setNames(split(fdf2, seq(nrow(fdf2))), rownames(fdf2))
  df.list <- lapply(df.list,FUN = function(x){
    text <- x[[1]]
    text <- strsplit(text,split = ",")
  })
  df.list <- lapply(df.list,unlist)
}


#'Computer for overlaping matrix based on given data
#'
#'@param pathwayid a vector that has all ramp pathway id in it
#'@param pathwaysWithAnalytes a list that has pathways ramp Id as name, analytes (ramp compound id only
#' or gene id only) as content.
#' @param methods must be in c('balanced','weighted') to determine which way to calculate this matrix
compute_overlap_matrix <- function(pathwayid,
                                   pathwaysWithAnalytes,
                                   methods){
  if(!(methods %in% c('balanced','weighted')))
    stop('Wrong option for the input')
  metabolite_result <- matrix(NA,nrow = length(pathwayid),ncol = length(pathwayid))
  colnames(metabolite_result) <- pathwayid
  rownames(metabolite_result) <- pathwayid
  # First method compute intersection over the union
  if(methods == 'balanced'){
    for(i in 1:length(pathwayid)){
      id <- pathwayid[i]
      cid <- pathwaysWithAnalytes[[i]]
      for (j in 1:length(pathwayid)) {
        if(is.na(metabolite_result[i,j])){
          if(i==j){
            metabolite_result[i,j] <- 1
          }else{
            cid2 <- pathwaysWithAnalytes[[j]]
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
  }else if (methods == 'weighted'){
    # second method
    for(i in 1:length(pathwayid)){
      id <- pathwayid[i]
      cid <- pathwaysWithAnalytes[[i]]
      for (j in 1:length(pathwayid)) {
        if(is.na(metabolite_result[i,j])){
          if(i==j){
            metabolite_result[i,j] <- 1
          }else{
            cid2 <- pathwaysWithAnalytes[[j]]
            shared_metabolite <- unique(intersect(cid,cid2))
            total <- unique(union(cid,cid2))
            metabolite_result[i,j] <- length(shared_metabolite)/length(unique(cid2))
            print(metabolite_result[i,j])
            if(is.na(metabolite_result[j,i])){
              metabolite_result[j,i] <- length(shared_metabolite)/length(unique(cid))
            }
          }
        }
        print(paste("Compute for ",i,",",j))
      }
    }
  }
  
  return(metabolite_result)
} # end for the functions


con <- dbConnect(MySQL(),
                 user = 'root',
                 dbname='Your-database-name',
                 password = 'Your-password',
                 host = 'localhost')

pathways<- dbGetQuery(con,'select * from pathway;')



dbname <- unique(pathways$type)

# pathwayInHmdb <- pathways[pathways$type == 'hmdb',]
pathwayInKegg <- pathways[pathways$type == 'kegg',]
pathwayInWiki <- pathways[pathways$type == 'wiki',]
pathwayInReac <- pathways[pathways$type == 'reactome',]

# define the minimum metabolites/genes
min_analyte <- 10

# Store Compound Ids in List
# listOfHmdbC <- findAnalyteHasPathway(pathwayInHmdb$pathwayRampId)
listOfKeggC <- findAnalyteHasPathway(pathwayInKegg$pathwayRampId,GC = 'both',n = min_analyte)
listOfWikiC <- findAnalyteHasPathway(pathwayInWiki$pathwayRampId,GC = 'both',n = min_analyte)
listOfReacC <- findAnalyteHasPathway(pathwayInReac$pathwayRampId,GC = 'both',n = min_analyte)

pathwayid <- c(#names(listOfHmdbC),
                names(listOfKeggC),
                names(listOfWikiC),
                names(listOfReacC))

pathToanalC <- do.call(c,list(#listOfHmdbC,
                              listOfKeggC,
                              listOfWikiC,
                              listOfReacC))
metabolite_result <- compute_overlap_matrix(pathwayid = pathwayid,
                                            pathwaysWithAnalytes =  pathToanalC,
                                            methods = 'balanced')
