args <- commandArgs(TRUE)
library(VennDiagram)
library(grDevices)

threeVenn <- function(wiki,
                      reactome,
                      kegg,
                      wr,
                      rk,
                      kr,
                      rwk){
    num1 = wiki
    num2 = reactome
    num3 = kegg
    num12 = wr
    num23 = rk
    num13 = kr
    num123 = rwk

    wikiPercent = round(wiki/(wiki+ reactome + kegg)*100, digits=0)
    reactomePercent = round(reactome/(reactome + kegg + wiki)*100, digits=0)
    keggPercent = round(kegg/((reactome + kegg + wiki))*100, digits=0)
    keggpercent = paste0("kegg (", keggPercent, "%)")
    wikipercent = paste0("wiki (", wikiPercent, "%)")
    reactomepercent = paste0("reactome (", reactomePercent, "%)")

    venn.plot <- draw.triple.venn(area1 = num1, 
                              area2 = num2, 
                              area3 = num3, 
                              n12 = num12, 
                              n23 = num23, 
                              n13 = num13, 
                              n123 = num123,
                              category = c(wikipercent, reactomepercent, keggpercent),
                              fill = c("orange", "red", "green"),
                              lty = "dashed",
                              cex = 2,
                              cat.cex = 2,
                              cat.col = c("orange", "red", "green"))
    
    name = paste0("../misc/output/", "threeVennApoptosis", ".pdf")
    pdf(file = name, width = 20, height = 20)
    grid.draw(venn.plot)
    dev.off()
}

threeVenn(as.numeric(args[1]),
         as.numeric(args[2]),
         as.numeric(args[3]),
         as.numeric(args[4]),
         as.numeric(args[5]),
         as.numeric(args[6]),
         as.numeric(args[7]))




