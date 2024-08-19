'''
Created on Nov 6, 2020

@author: braistedjc
'''

class NamedSet:
    def __init__(self, name):
        self.inner_set = set()
        self.name = name

    def add(self, str):
        self.inner_set.add(str)

    def update(self, other_set):
        self.inner_set.update(other_set.inner_set)

    def __hash__(self):
        return id(self.inner_set)

    def __eq__(self, other):
        return id(self.inner_set) == id(other.inner_set)

    def __repr__(self):
        return f"NamedSet(name={self.name}, inner_set={self.inner_set})"

class MetaboliteList(object):
    '''
    Container list class holding metabolites. The list class supports access and set methods as well as
    utility methods to extract metabolite records.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.metaboliteSourceIdList = dict()

        self.inchikeyPrefixToMetab = dict()

        self.sourceSummary = dict()


    def contains(self, metabolite):
        """Utility method returning true if the metabolite exists in the list or false if
        the metabolite does not exist in the list based on source id.
        """
        return (self.metaboliteSourceIdList[metabolite.sourceId] != None)


    def getMetaboliteBySourceId(self, sourceId):
        """
        returns a metabolite based on the source id
        """
        return self.metaboliteSourceIdList.get(sourceId)


    def addMetabolite(self, metabolite):
        """
        Adds a metabolite to the list
        """
        self.metaboliteSourceIdList[metabolite.sourceId] = metabolite


    def addMetaboliteByAltId(self, id, metabolite):
        """
        Adds a metabolite connected to a supplied id.
        This permits multiple ids to point to a shared metabolite.
        """
        self.metaboliteSourceIdList[id] = metabolite


    def length(self):
        """
        Returns the number of source ids, not distinct metabolites.
        """
        return len(self.metaboliteSourceIdList)

    def getUniqueMetabolites(self):
        """
        Returns the set of unique metabolite entities as a list.
        """
        return list(set(self.metaboliteSourceIdList.values()))

    def getAllMetabolites(self):
        return list(self.metaboliteSourceIdList.values())

    def generateMetaboliteSourceStats(self, sourceList):
        """
        Returns a map of source to count of metabolites from each source.
        """
        for source in sourceList:
            self.sourceSummary[source.sourceName] = 0

        mets = self.getUniqueMetabolites()

        for met in mets:
            for source in met.sources:
                self.sourceSummary[source] = self.sourceSummary[source] + 1

        return self.sourceSummary


    def printChemPropSummaryStats(self):
        """
        Utility to print chemical property summary statistics for chemical properties.
        """
        mets = self.getUniqueMetabolites()

        haveMolCount = 0
        molRecords = 0

        for met in mets:
            if len(met.chemPropsMolecules) > 0:
                haveMolCount = haveMolCount + 1

                for source in met.chemPropsMolecules:
                    molRecords = molRecords + len(met.chemPropsMolecules[source])

        print("\nChemistry Property Stats")
        print("Tot Mets: " + str(len(mets)))
        print("Tot Mets with ChemProps: " + str(haveMolCount))
        print("Tot Molecules records: " + str(molRecords))



    def get_inchi_keys(self, metabolite):
        median_mw = metabolite.get_median_mw()
        if median_mw >= 500:
            met1_inchis = metabolite.getInchiPrefixes()
        else:
            met1_inchis = metabolite.getInchiKeyDuplexes()
        return met1_inchis

    def collapseMetsOnInchiKeyPrefix(self):
        metabolites = self.getUniqueMetabolites()
        inchi_set_index = dict()

        for met in metabolites:
            keys = self.get_inchi_keys(met)

            current_set = NamedSet(met.rampId)

            for key in keys:
                current_set.add(key)
                if key in inchi_set_index:
                    existing_set = inchi_set_index[key]

                    existing_set.update(current_set)
                    for merged_key in current_set.inner_set:
                        inchi_set_index[merged_key] = existing_set

                    current_set = existing_set
                else:
                    inchi_set_index[key] = current_set

        for met in metabolites:
            keys = self.get_inchi_keys(met)
            if len(keys) > 0:
                inchi_set = inchi_set_index[next(iter(keys))]
                updated_id = inchi_set.name
                if met.rampId != updated_id:
                    met.rampId = updated_id

