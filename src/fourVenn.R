args <- commandArgs(TRUE)
library(VennDiagram)
library(grDevices)

fourVenn <- function(kegg, 
                     hmdb, 
                     reactome, 
                     wiki, 
                     kh,
                     kr,
                     kw,
                     hr,
                     hw,
                     rw,
                     khr,
                     krw,
                     hrw,
                     khw,
                     khrw,
                     type) {

    num1 = kegg
    num2 = hmdb
    num3 = reactome
    num4 = wiki
    num12 = kh
    num13 = kr
    num14 = kw
    num23 = hr
    num24 = hw
    num34 = rw
    num123 = khr
    num134 = krw
    num234 = hrw
    num124 = khw
    num1234 = khrw
    
    khrPercent = khr - khrw
    krwPercent = krw - khrw
    hrwPercent = hrw - khrw
    khwPercent = khw - khrw
    kwPercent = kw - krwPercent - khwPercent - khrw
    hwPercent = hw  - hrwPercent - khwPercent - khrw
    rwPercent = rw - krwPercent - hrwPercent - khrw
    khPercent = kh - khrPercent - khwPercent - khrw
    hrPercent = hr - khrPercent - hrwPercent - khrw
    krPercent = kr - khrPercent - krwPercent - khrw
    wikiPercent = wiki - kwPercent - hwPercent - rwPercent - krwPercent - hrwPercent - khwPercent -  khrw
    hmdbPercent = hmdb - khPercent - hrPercent - hwPercent - khrPercent - hrwPercent - khwPercent - khrw
    reactomePercent = reactome - krPercent - hrPercent - rwPercent - khrPercent - krwPercent - hrwPercent - khrw
    keggPercent = kegg - khPercent - krPercent - kwPercent - krwPercent - khwPercent - khrPercent - khrw 
    
    keggpercent = round(((keggPercent)/kegg)*100, digits=1)
    hmdbpercent = round(((hmdbPercent)/hmdb)*100, digits=1)
    reactomepercent = round((reactomePercent)/reactome*100, digits=1)
    wikipercent = round(((wikiPercent)/wiki)*100, digits=1)
    keggpercent = paste0("kegg (", keggpercent, "%)")
    hmdbpercent = paste0("hmdb (", hmdbpercent, "%)")
    reactomepercent = paste0("reactome (", reactomepercent, "%)")
    wikipercent = paste0("wiki (", wikipercent, "%)")
    
    venn.plot <- draw.quad.venn(
        area1 = num1,
        area2 = num2,
        area3 = num3,
        area4 = num4,
        n12 = num12,
        n13 = num13,
        n14 = num14,
        n23 = num23,
        n24 = num24,
        n34 = num34,
        n123 = num123,
        n124 = num124,
        n134 = num134,
        n234 = num234,
        n1234 = num1234,
        category = c(keggpercent, hmdbpercent, reactomepercent, wikipercent),
        fill = c("orange", "red", "green", "blue"),
        lty = "dashed",
        cex = 2,
        cat.cex = 2,
        cat.col = c("orange", "red", "green", "blue"))
        
        name = paste0("../misc/output/", "fourVenn", type, ".pdf")
        pdf(file = name, width = 20, height = 20)
        grid.draw(venn.plot)
        dev.off()
}

fourVenn(as.numeric(args[1]),
         as.numeric(args[2]),
         as.numeric(args[3]),
         as.numeric(args[4]),
         as.numeric(args[5]),
         as.numeric(args[6]),
         as.numeric(args[7]),
         as.numeric(args[8]),
         as.numeric(args[9]),
         as.numeric(args[10]),
         as.numeric(args[11]),
         as.numeric(args[12]),
         as.numeric(args[13]),
         as.numeric(args[14]),
         as.numeric(args[15]),
         args[16])

