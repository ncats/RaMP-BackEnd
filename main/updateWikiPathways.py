from updateSQL import RampUpdater
import pickle as pk
def getUpdateObjectPkl():
    pkf = open('../misc/output/wikipathwayRdfPk.pkl','rb')
    wp = pk.load(pkf)
    print(wp.pathwayWithMetabolitesDictionary)
    pkf.close()
    rp = RampUpdater(wp)
    rp.checkNewAnalyteEntry('compound')
    rp.checkNewPathwayEntry()
    pkf2 = open('../misc/output/updateObject328.pkl','wb')
    pk.dump(rp,pkf2,pk.HIGHEST_PROTOCOL)
    pkf2.close()
    print('New pathway: {}'.format(len(rp.newRampPathway)))
    
def updatingRamp():
    pkf = open('../misc/output/updateObject328.pkl','rb')
    rp = pk.load(pkf)
    pkf.close()
    rp2 = RampUpdater(rp)
    rp2.newRampCompound = rp.newRampCompound
    rp2.newRampPathway = rp.newRampPathway
    rp2.writeToRamp('wiki')
if __name__ == '__main__':
    updatingRamp()