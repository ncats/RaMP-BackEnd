import subprocess
class getStatistics():
    
    '''
    The purpose the statistics class is to get an overall "feel" of the information contained in the RaMP database.
    The word "statistics" refers to the things like the number of genes in the database, the number of pathways, etc. 
    
    '''
    
    
    def analyteOverlaps(self, rampIdInWhichDatabases, rampIdDictionary, analyteType, writeToFiles = False):
        
        '''
        The purpose of the function analyteOverlaps is to compare information among the four databases that make 
        up the RaMP databases -- wikipathways, reactome, kegg, hmdb. There may be overlaps between two ormore 
        of the databases. For example, one metabolite may be present in two databaase. The metabolite in  database 1
        may be identified by the same metabolite id in both databases, or they may be different. Differences in IDs are
        resolved in the ID conversion class -- by the time this function is run all known ids for a metabolite are contained 
        in the rampIdDictionary dictionary, which links the RAMPID to every metabolite id.
        
        The output of the function is a four-way venn diagram showing analyte overlaps.
        
        Although the example above refers only to metabolites, this function works for either genes or metabolites. 
        However, they must be submitted separately. 
        
        Analyte means gene OR metabolite. The function works on either genes OR metabolites. 
        
        param dict rampIdInWhichDatabases: dictionary where the rampID is the key and a set of databases is the value. This 
        dictionary keeps track of which databases a RAMPID can be found in. For example, hmdb and kegg. This dictionary is created 
        in the writeToSQL class.
        param dict rampIdDictionary: dictionary containing either metabolites or genes, with the analyte databaseID as the key  
        and the analyte RAMPID as the value. This dictionary is created in the writeToSQL class. 
        param str analyteType: "Gene" or "Compound"
        param bool writeToFiles: A file can be created where each line is an analyte. The first column is a list of identifies
        for analyte and the second column is the databases the analyte is found in (abbreviated, k for kegg, r for reactome, etc).
        This parameter is default off since this is time-consuming.
        '''
        
        statisticsOutFile = open("../misc/output/statisticsAnalyteOverlaps" + str(analyteType) + ".txt", 'wb')
        
        
        kegg = 0
        hmdb = 0
        reactome = 0
        wiki = 0
        kh = 0
        kr = 0
        kw = 0
        hr = 0
        hw = 0
        rw = 0
        khr = 0
        krw = 0
        hrw = 0
        khw = 0
        khrw = 0
        
        keggSet = set(["kegg"])
        hmdbSet = set(["hmdb"])
        reactomeSet = set(["reactome"])
        wikiSet = set(["wiki"])
        khSet = set(["kegg", "hmdb"])
        krSet = set(["kegg", "reactome"])
        kwSet = set(["kegg", "wiki"])
        hrSet = set(["hmdb", "reactome"])
        hwSet = set(["hmdb", "wiki"])
        rwSet = set(["reactome", "wiki"])
        khrSet = set(["kegg", "hmdb", "reactome"])
        krwSet = set(["kegg", "reactome", "wiki"])
        hrwSet = set(["hmdb", "reactome", "wiki"])
        khwSet = set(["kegg", "hmdb", "wiki"])
        
        khrwSet = set(["kegg", "hmdb", "wiki", "reactome"])
        
        
        if writeToFiles:
            rampIDtoanalyteIdDictionary = {}
            num = 0
            for rampID in rampIdInWhichDatabases:
                print("processing analyte number: " + str(num) + "/" + str(len(rampIdInWhichDatabases)))
                num = num + 1
                analyteList = []
                for analyteID in rampIdDictionary:
                    if rampIdDictionary[analyteID] == rampID:
                        analyteList.append(analyteID)
            
                rampIDtoanalyteIdDictionary[rampID] = analyteList
        
        for rampID in rampIdInWhichDatabases:
            
            rampSet = rampIdInWhichDatabases[rampID]
            
            if rampSet == keggSet:
                kegg = kegg + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",kegg\n")
            
            elif rampSet == hmdbSet:
                hmdb = hmdb + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",hmdb\n")
            
            elif rampSet == reactomeSet:
                reactome = reactome + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",reactome\n")
            
            elif rampSet == wikiSet:
                wiki = wiki + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",wiki\n")
            
            elif rampSet == khSet:
                kh = kh + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",kh\n")
            
            elif rampSet == krSet:
                kr = kr + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",kr\n")
            
            elif rampSet == kwSet:
                kw = kw + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",kw\n")
            
            elif rampSet == hrSet:
                hr = hr + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",hr\n")
            
            elif rampSet == hwSet:
                hw = hw + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",hw\n")
            
            elif rampSet == rwSet:
                rw = rw + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",rw\n")
            
            elif rampSet == khrSet:
                khr = khr + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",khr\n")
            
            elif rampSet == krwSet:
                krw = krw + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",krw\n")
            
            elif rampSet == hrwSet:
                hrw = hrw + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",hrw\n")
            
            elif rampSet == khwSet:
                khw = khw + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",khw\n")
            
            elif rampSet == khrwSet:
                khrw = khrw + 1
                
                if writeToFiles:
                    value = rampIDtoanalyteIdDictionary[rampID]
                    statisticsOutFile.write(":".join(value).encode("utf-8") + b",khrw\n")
        
        
        if writeToFiles:
            statisticsOutFile.close()
        
        wiki = wiki + kw + hw + rw + krw + hrw + khw +  khrw
        hmdb = hmdb + kh + hr + hw + khr + hrw + khw + khrw
        reactome = reactome + kr + hr + rw + khr + krw + hrw + khrw
        kegg = kegg + kh + kr + kw + krw + khw + khr + khrw 
        
        kw = kw + krw + khw +  khrw
        hw = hw + hrw + khw +  khrw
        rw = rw + krw + hrw +  khrw
        kh = kh + khr + khw + khrw
        hr = hr + khr + hrw + khrw
        kr = kr + khr + krw + khrw
        
        krw = krw +  khrw
        hrw = hrw +  khrw
        khw = khw +  khrw
        khr = khr + khrw
        
        print("kegg: " + str(kegg))
        print("hmdb: " + str(hmdb))
        print("reactome: " + str(reactome))
        print("wiki: " + str(wiki))
        print("kh: " + str(kh))
        print("kr: " + str(kr))
        print("kw: " + str(kw))
        print("hr: " + str(hr))
        print("hw: " + str(hw))
        print("rw: " + str(rw))
        print("khr: " + str(khr))
        print("krw: " + str(krw))
        print("hrw: " + str(hrw))
        print("khw: " + str(khw))
        print("khrw: " + str(khrw))
        
        bash = ""
        
        for each in [kegg, hmdb, reactome, wiki, kh, kr, kw, hr, hw, rw, khr, krw, hrw, khw, khrw]:
           bash = str(bash) + str(each) + " " 
            
        print(bash)
        
        bash = "Rscript ../__init__/fourVenn.R " + str(bash) + str(analyteType)
        
        print(bash)
        subprocess.call(bash, shell=True)
        
    def databaseContent(self, pathwayDictionary, 
                     pathwayCategory,
                     metabolitesWithPathwaysDictionary,
                     metabolitesWithSynonymsDictionary,
                     metaboliteIDDictionary,
                     pathwaysWithGenesDictionary,
                     geneInfoDictionary,
                     biofluidLocation,
                     biofluid,
                     cellularLocation,
                     cellular,
                     pathwayOntology,
                     endoExoDictionary,
                     database):
        
        '''
        The purpose of databaseContent is to get overall metrics for the individual databases that make up RaMP. 
        These metrics will be printed out to the screen. They are:
        
        1. Number of genes
        2. Number of metabolites
        3. Number of pathways
        4. Biofluids (if applicable, only found in hmdb)
        5. Cellular location (if application, only found in hmdb)
        
        param dict pathwayDictionary: see class for database 
        param dict pathwayCategory: see class for database
        param dict metabolitesWithPathwaysDictionary: see class for database
        param dict metabolitesWithSynonymsDictionary: see class for database
        param dict metaboliteIDDictionary: see class for database
        param dict pathwaysWithGenesDictionary: see class for database
        param dict geneInfoDictionary: see class for database
        param dict biofluidLocation: see class for database (may be empty for some classes)
        param dict biofluid: see class for database (may be empty for some classes)
        param dict cellularLocation: see class for database (may be empty for some classes)
        param dict cellular: see class for database (may be empty for some classes)
        param dict pathwayOntology: see class for database (may be empty for some classes)
        param dict endoExoDictionary: see class for database (may be empty for some classes)
        param str  database: name of the database (e.g. "kegg") 
        '''
            
        numPathways = len(pathwayDictionary)
        numMetabolites = len(metaboliteIDDictionary)
        numGenes = len(geneInfoDictionary)
       
        print(database)
        print("Number of pathways: " + str(numPathways))
        print("Number of metabolites: " + str(numMetabolites))
        print("Number of genes: " + str(numGenes))
       
       
        if len(endoExoDictionary) > 0:
            food = 0
            endogenous = 0
            drugmetabolite = 0
            toxinpollutant = 0 
            drug = 0
            microbial = 0
            plant = 0
            cosmetic = 0
            drugorsteroid = 0
            exogenous = 0
            micorbial = 0  
            
            for key in endoExoDictionary:
                value = endoExoDictionary[key]
                for each in value:
                    if each == "Food":
                        food += 1
                    if each == "Endogenous":
                        endogenous += 1
                    if each == "Drug metabolite":
                        drugmetabolite += 1
                    if each == "Toxin/Pollutant":
                        toxinpollutant += 1
                    if each == "Drug":
                        drug += 1
                    if each == "Microbial":
                        microbial += 1
                    if each == "Plant":
                        plant += 1
                    if each == "Cosmetic":
                        cosmetic += 1
                    if each == "Drug or steroid metabolite":
                        drugorsteroid += 1
                    if each == "Exogenous":
                        exogenous += 1
                    if each == "Micorbial":
                        micorbial += 1
                
            
            print("***Exo/Endo***")
            print("Food: " + str(food))
            print("Endogenous: " + str(endogenous))
            print("Drug metabolite: " + str(drugmetabolite))
            print("Toxin/Pollutant: " + str(toxinpollutant))
            print("Drug: " + str(drug))
            print("Microbial: " + str(microbial))
            print("Plant: " + str(plant))
            print("Cosmetic: " + str(cosmetic))
            print("Drug or steroid metabolite: " + str(drugorsteroid))
            print("Exogenous: " + str(exogenous))
            print("Micorbial: " + str(micorbial))
        
        if len(biofluidLocation) > 0: 
            feces = 0
            amniotic = 0
            breastmilk = 0
            urine = 0
            blood = 0
            aqhumor = 0
            sweat = 0
            prostate = 0
            pericard = 0
            saliva = 0
            semen = 0
            tears = 0
            bile = 0
            ascites = 0
            cerebrospinal = 0
            cytoplasm = 0
            lymph = 0 
            
            for key in biofluidLocation:
                value = biofluidLocation[key]
                for each in value:
                    if each == "Feces":
                        feces += 1
                    if each == "Amniotic Fluid":
                        amniotic += 1
                    if each == "Breast Milk":
                        breastmilk += 1
                    if each == "Urine":
                        urine += 1
                    if each == "Blood":
                        blood += 1
                    if each == "Aqueous Humour":
                        aqhumor += 1
                    if each == "Sweat":
                        sweat += 1
                    if each == "Prostate Tissue":
                        prostate += 1
                    if each == "Pericardial Effusion":
                        pericard += 1
                    if each == "Saliva":
                        saliva += 1
                    if each == "Semen":
                        semen += 1
                    if each == "Tears":
                        tears += 1
                    if each == "Bile":
                        bile += 1
                    if each == "Ascites Fluid":
                        ascites += 1
                    if each == "Cerebrospinal Fluid (CSF)":
                        cerebrospinal += 1
                    if each == "Cellular Cytoplasm":
                        cytoplasm += 1
                    if each == "Lymph":
                        lymph += 1
        
            print("***biofluids***")
            print("Feces: " + str(feces))
            print("Amniotic Fluid: " + str(amniotic))
            print("Breast Milk: " + str(breastmilk))
            print("Urine: " + str(urine))
            print("Blood: "+ str(blood))
            print("Aqueous Humour: " + str(aqhumor))
            print("Sweat: " + str(sweat))
            print("Prostate Tissue: " + str(prostate))
            print("Pericardial Effusion: " + str(pericard))
            print("Saliva: " + str(saliva))
            print("Semen: " + str(semen))
            print("Tears: " + str(tears))
            print("Bile: " + str(bile))
            print("Ascites Fluid: " + str(ascites))
            print("Cerebrospinal Fluid (CSF): " + str(cerebrospinal))
            print("Cellular Cytoplasm: " + str(cytoplasm))
            print("Lymph: " + str(lymph))
        
        
        if len(cellularLocation) > 0:    
            cytologP = 0
            mito = 0
            extracell = 0
            nuc = 0
            endoretic = 0
            innermito = 0
            cyto = 0
            lyso = 0
            perox = 0
            membraneloP = 0
            golgi = 0
            membrane = 0
            microsomes = 0
            
    
            
            for key in cellularLocation:
                value = cellularLocation[key]
                for each in value:
                    if each == "Cytoplasm (predicted from logP)":
                        cytologP += 1
                    if each == "Mitochondria":
                        mito += 1
                    if each == "Extracellular":
                        extracell += 1
                    if each == "Nucleus":
                        nuc += 1
                    if each == "Endoplasmic reticulum":
                        endoretic += 1
                    if each == "Inner mitochondrial membrane":
                        innermito += 1
                    if each == "Cytoplasm":
                        cyto += 1
                    if each == "Lysosome":
                        lyso += 1
                    if each == "Peroxisome":
                        perox += 1
                    if each == 'Membrane (predicted from logP)':
                        membraneloP += 1
                    if each == 'Golgi apparatus':
                        golgi += 1
                    if each == 'Membrane':
                        membrane += 1
                    if each == 'Microsomes':
                        microsomes += 1
                        
    
           
            
            
            print("***cellular***")
            print("Cytoplasm (predicted from logP): " + str(cytologP))
            print("Mitochondria: " + str(mito))
            print("Extracellular: " + str(extracell))
            print("Nucleus: " + str(nuc))
            print("Endoplasmic reticulum: " + str(endoretic))
            print("Inner mitochondrial membrane: " + str(innermito))
            print("Cytoplasm: " + str(cyto))
            print("Lysosome: " + str(lyso))
            print("Peroxisome: " + str(perox))
            print("Membrane (predicted from logP): " + str(membraneloP))
            print("Golgi apparatus: " + str(golgi))
            print("Membrane: " + str(membrane))
            print("Microsomes: " + str(microsomes))
            
            
            
    def Apoptosis(self, rampIdDictionary, 
                  pathwaysWithGenesDictionaryWiki,
                  pathwaysWithGenesDictionaryKegg,
                  pathwaysWithGenesDictionaryReactome):
        
                '''
                The purpose of the apoptosis function is to identify the genes present in the pathway called "Apoptosis"
                in 3 of the 4 databases in ramp: "wikipathways", "reactome", and "hmdb" and see if there is overlap 
                among the genes. This overlap will be represent in a 3-way venn diagram. 
                
                This figure is a replication of the same figure done in a published paper and is a "quality check" to 
                ensure that the various steps involved in creation of the RaMP database such as getting the data and id 
                conversion have proceeded smoothly. 
                
                The paper is here: https://www.nature.com/nprot/journal/v11/n10/full/nprot.2016.117.html
                
                param dict rampIdDictionary: dictionary containing genes, with the gene databaseID as the key  
        and the analyte RAMPID as the value. This dictionary is created in the writeToSQL class.
                param dict pathwaysWithGenesDictionaryWiki: pathways in wikipathways
                param dict pathwaysWithGenesDictionaryKegg: pathways in kegg
                param dict pathwaysWithGenesDictionaryReactome: pathways in reactome
                
                
                '''
        
                wikiGenes = pathwaysWithGenesDictionaryWiki["WP254"]
                reactomeGenes = pathwaysWithGenesDictionaryReactome["R-HSA-109581"]
                keggGenes = pathwaysWithGenesDictionaryKegg["04210"]
                
                wikiGenesRamp = set()
                reactomeGenesRamp = set()
                keggGenesRamp = set()
                
                for item in wikiGenes:
                    rampID = rampIdDictionary[item]
                    wikiGenesRamp.add(rampID)
                
                for item in reactomeGenes:
                    rampID = rampIdDictionary[item]
                    reactomeGenesRamp.add(rampID)
                
                for item in keggGenes:
                    rampID = rampIdDictionary[item]
                    keggGenesRamp.add(rampID)
                
                
                
                w=len(wikiGenesRamp)
                r=len(reactomeGenesRamp)
                k=len(keggGenesRamp)
                wr=len(wikiGenesRamp.intersection(reactomeGenesRamp))
                rk=len(reactomeGenesRamp.intersection(keggGenesRamp))
                kw=len(keggGenesRamp.intersection(wikiGenesRamp))
                rwk=len(keggGenesRamp.intersection(wikiGenesRamp).intersection(reactomeGenesRamp))
                

                bash = ""
        
                for each in [w, r, k, wr, rk, kw, rwk]:
                    bash = str(bash) + str(each) + " " 
            
                print(bash)
        
                bash = "Rscript ../__init__/threeVenn.R " + str(bash)
        
                print(bash)
                subprocess.call(bash, shell=True)
                
                
        
            
            
            