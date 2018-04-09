from keggBacterial import keggBacterial

if __name__ == '__main__':
    kegg = keggBacterial()
    kegg.getAllSpecies()
    kegg.getSpeciesPathwayFile()
    kegg.getAllPathwaysFromMicrobial()
    kegg.getAllLinkedCpdFromPathways()
    kegg.getAllCpdFromMapFile()