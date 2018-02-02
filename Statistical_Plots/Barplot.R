library(RMySQL)

# Define the functions
# find the analytes' ramp ids that have pathway associated 
finderPathways <- function(SourceId,rampIdHasPathway){
  hmdbSourceId <- unique(SourceId$V1)
  hmdbSourceId <- sapply(hmdbSourceId,shQuote)
  hmdbSourceId <- paste(hmdbSourceId,collapse = ",")
  query <- paste0("select * from source where sourceId in (",hmdbSourceId,");")
  df <- dbGetQuery(con,query)
  #return(df)
  rampId <- intersect(df$rampId,rampIdHasPathway)
}
# Get the analytes-pathway relations from database
finderPathways2 <- function(Ids){
  Ids <- unique(Ids)
  Ids <- sapply(Ids,shQuote)
  Ids <- paste(Ids,collapse = ",")
  query <- paste0("select * from analytehasPathway where rampId in (",Ids,");")
  df <- dbGetQuery(con,query)
}
# return the dataframe that has database name as column, ramp metabolite id as row
# each cell is the number of metabolites that have this pathway.
finderPathways3 <- function(dfOfIds){
  cid <- unique(dfOfIds$rampId)
  pid <- unique(dfOfIds$pathwayRampId)
  pathSource <- c("hmdb","kegg","wiki","reactome")
  output <- data.frame(metabolites = NULL,totpathways = NULL,hmdb=NULL,
                       kegg =NULL,wiki=NULL,reac = NULL) 
  for(id in cid){
    print(id)
    pathways <- dfOfIds[dfOfIds$rampId == id,]
    numOfPathway <- vector()
    for(db in pathSource){
      numOfPathway <- c(numOfPathway,nrow(unique(pathways[pathways$pathwaySource == db,])))
      # print(pathways[pathways$type == db,])
      # Sys.sleep(3)
      
    }
    output <- rbind(output,data.frame(
      metabolites = id,totpathways = nrow(pathways),
      hmdb = numOfPathway[1],
      kegg = numOfPathway[2],
      wiki = numOfPathway[3],
      reac = numOfPathway[4]
    ))
  }
  return(output)
}
# End of functions

con <- dbConnect(MySQL(),
                 dbname ="mathelabramp",
                 host ="localhost",
                 password = "Ehe131224",
                 username = "root")
# Read all source id from Python output
# The data is output from Python code.
# Make sure the working directory has these files
hmdb <- read.table('data/hmdbmetaboliteIDDictionary.txt')
kegg <- read.table("data/keggmetaboliteIDDictionary.txt")
wiki <- read.table("data/wikimetaboliteIDDictionary.txt")
reac <- read.table("data/reactomemetaboliteIDDictionary.txt")

# Use unique ID in the first column to search
hmdb <- data.frame(V1 = as.character(unique(hmdb$V1)))
kegg <- data.frame(V1 = as.character(unique(kegg$V1)))
wiki <- data.frame(V1 = as.character(unique(wiki$V1)))
reac <- data.frame(V1 = as.character(unique(reac$V1)))

query <- "select * from analytehaspathway;"
analyteHasPathway <- unique(dbGetQuery(con,query))


rampIdHasPathway <- unique(analyteHasPathway$rampId)
# Find all unique metabolites from each database
rampIdForHMDB <- finderPathways(hmdb,rampIdHasPathway)
rampIdForKEGG <- finderPathways(kegg,rampIdHasPathway)
rampIdForWIKI <- finderPathways(wiki,rampIdHasPathway)
rampIdForREAC <- finderPathways(reac,rampIdHasPathway)
# Find how many metabolites for each metabolites

pathwaysForHMDB <- finderPathways2(rampIdForHMDB)
resultHMDB <- finderPathways3(pathwaysForHMDB)

pathwaysForKEGG <- finderPathways2(rampIdForKEGG)
resultKEGG <- finderPathways3(pathwaysForKEGG)

pathwaysForWIKI <- finderPathways2(rampIdForWIKI)
resultWIKI <- finderPathways3(pathwaysForWIKI)

pathwaysForREAC <- finderPathways2(rampIdForREAC)
resultREAC <- finderPathways3(pathwaysForREAC)

# Generate bar plot.
library(plotly)
# Initilize plot data for each database
plot.hmdb <- resultHMDB[order(resultHMDB$hmdb,decreasing = T) & resultHMDB$hmdb!=0,]
plot.kegg <- resultKEGG[order(resultKEGG$kegg,decreasing = T) & resultKEGG$kegg!=0,]
plot.wiki <- resultWIKI[order(resultWIKI$wiki,decreasing = T) & resultWIKI$wiki!=0,]
plot.reac <- resultREAC[order(resultREAC$reac,decreasing = T) & resultREAC$reac!=0,]
# Store them in a list for loop
plot.data <- list()
plot.data[[1]] <- plot.hmdb
plot.data[[2]] <- plot.kegg
plot.data[[3]] <- plot.wiki
plot.data[[4]] <- plot.reac

# Find synonyms for each database 
database <- c("hmdb",'kegg','wiki','reac')
#
for(i in 1:4){
  rampId <- as.character(unique(plot.data[[i]]$metabolites))
  rampId <- sapply(rampId,shQuote)
  rampId <- paste(rampId,collapse = ",")
  query <- paste("select * from analytesynonym where rampId in(",rampId,") group by rampId;")
  synonymdf1 <- dbGetQuery(con,query)
  plot.data[[i]] <- merge(plot.data[[i]],synonymdf1[,1:2],by.x = "metabolites",by.y="rampId",all.x = T)
}
for(i in 1:4){
  plot.data[[i]] <- plot.data[[i]][order(plot.data[[i]][,i+2],decreasing = T),]
}


# Set up height,width, and margins for plot
m <- list(
  l = 100,
  r = 100,
  t = 100,
  b = 100
)
h <- 1000
w <- 1000

p_hmdb <- plot_ly(x = 1:length(plot.data[[1]]$Synonym),
                  y = plot.data[[1]]$hmdb,
                  type = "bar",
                  height = h,
                  width = w) %>%
  layout(title = "HMDB metabolites vs. Pathways",
         margin = m,
         font = list(
           size = 24
         ),
         yaxis = list(
           title = "Number of pathways metabolites are involved in",
           font = list(
             size = 24
           )
         ),
         xaxis = list(
           title = paste("Total", nrow(resultHMDB),"unique metabolites have pathways from HMDB"),
           showticklabels = FALSE,
           font = list(
             size = 24
           )
         ))
# For loop add annotation to the plot.
# ax,ay is the cordinates of the annotation, they're hard coded
# to find good position 
for(i in 1:5){
  df <- plot.data[[1]]
  ay <- 40
  if(i == 4){
    ay <- -15
    ax <- 40
  }else if(i ==5){
    ay <- 15
    ax <- 80
  } else{
    ay <- 40
    ax <- 40
  }
  p_hmdb <- p_hmdb %>%
    add_annotations(
      text = df$Synonym[i],
      x = i,
      y = df$hmdb[i],
      ax = 40,
      ay = ay,
      xanchor = "left",
      yanchor = "bottom"
    )
}

p_hmdb


p_kegg <- plot_ly(x = 1:length(plot.data[[2]]$Synonym),
                  y = plot.data[[2]]$kegg,
                  type = "bar",
                  height = h,
                  width = w) %>%
  layout(title = "KEGG metabolites vs. Pathways",
         margin = m,
         font = list(
           size = 24
         ),
         yaxis = list(
           title = "Numbers of pathways metabolites are involved in",
           font = list(
             size = 24
           )
         ),
         xaxis = list(
           title = paste("Total", nrow(resultKEGG),"unique metabolites have pathways from KEGG"),
           showticklabels = FALSE,
           font = list(
             size = 24
           )
         )) 

for(i in 1:5){
  df <- plot.data[[2]]
  ax <- 20
  ay <- 20
  if(i == 2){
    ax <- 20
    ay <- -20
  } else{
    ax <- 20
    ay <- 20
  }
  p_kegg <- p_kegg %>%
    add_annotations(
      text = df$Synonym[i],
      x = i,
      y = df$kegg[i],
      ax = ax,
      ay = ay,
      xanchor = "left",
      yanchor = "bottom"
    )
}
p_kegg


p_wiki <- plot_ly(x = 1:length(plot.data[[3]]$Synonym),
                  y = plot.data[[3]]$wiki,
                  type = "bar",
                  height = h,
                  width = w) %>%
  layout(title = "Wikipathways metabolites vs. Pathways",
         margin = m,
         font = list(
           size = 24
         ),
         yaxis = list(
           title = "Numbers of Wiki pathways that metabolites are involved in",
           font = list(
             size = 24
           )
         ),
         xaxis = list(
           title = paste("Total", nrow(resultWIKI),"unique metabolites have pathways from WP"),
           showticklabels = FALSE,
           font = list(
             size = 24
           )
         )) 
for(i in 1:5){
  df <- plot.data[[3]]
  ax <- 20
  ay <- 20
  if(i == 2 ) {
    ax <- 20
    ay <- -20
  } else if (i == 3){
    ax <- 20
    ay <- 20
  } else if (i == 4) {
    ax <- 40
    ay <- 50
  } else if (i == 5) {
    ax <- 60
    ay <- 80
  }
  p_wiki <- p_wiki %>%
    add_annotations(
      text = df$Synonym[i],
      x = i,
      y = df$wiki[i],
      ax = ax,
      ay = ay,
      xanchor = "left",
      yanchor = "bottom"
    )
}


p_wiki


p_reac <- plot_ly(x = 1:length(plot.data[[4]]$reac),
                  y = plot.data[[4]]$reac,
                  type = "bar",
                  height = h,
                  width = w) %>%
  layout(title = "Reactome metabolites vs. Pathways",
         margin = m,
         font = list(
           size = 24
         ),
         yaxis = list(
           title = "Numbers of pathways metabolites are involved in",
           font = list(
             size = 24
           )
         ),
         xaxis = list(
           title = paste("Total", nrow(resultREAC),"unique metabolites have pathways from Reactome"),
           showticklabels = FALSE,
           font = list(
             size = 24
           )
         ))
for(i in 1:5){
  df <- plot.data[[4]]
  p_reac <- p_reac %>%
    add_annotations(
      text = df$Synonym[i],
      x = i,
      y = df$reac[i],
      ax = 20,
      ay = 20,
      xanchor = "left",
      yanchor = "bottom"
    )
}

p_reac