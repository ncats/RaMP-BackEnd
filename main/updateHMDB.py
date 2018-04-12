from updateSQL import RampUpdater,QualityControl
import pickle as pk
import time
from hmdbData import hmdbData
def getHmdbPkl():
    hmdb = hmdbData()
    hmdb.getEverything(True)
    hmdb_pkl = open('../misc/output/hmdbPkl.pkl','wb')
    pk.dump(hmdb,hmdb_pkl,pk.HIGHEST_PROTOCOL)
    hmdb_pkl.close()
    del hmdb_pkl
def getUpdateObjectPkl():
    pkf = open('../misc/output/hmdbPkl.pkl','rb')
    hmdb = pk.load(pkf)
    print(hmdb.__dir__())
    pkf.close()
    rp = RampUpdater(hmdb)
    rp.checkNewAnalyteEntry('compound')
    rp.checkNewAnalyteEntry('gene')
    rp.checkNewPathwayEntry()
    pkf2 = open('../misc/output/updateHMDBObject411.pkl','wb')
    pk.dump(rp,pkf2,pk.HIGHEST_PROTOCOL)
    pkf2.close()
def updatingRamp(db):
    pkf = open('../misc/output/updateObject328.pkl','rb')
    rp = pk.load(pkf)
    print('Unique metabolites: {} \nUnique genes: {}'.format(len(set(rp.newRampCompound.values())),
                                                              len(set(rp.newRampGene.values()))))
    
    time.sleep(1)
    pkf.close()
    rp2 = RampUpdater(rp)
    rp2.newRampCompound = rp.newRampCompound
    rp2.newRampPathway = rp.newRampPathway
    rp2.newRampGene = rp.newRampGene
    rp2.writeToRamp(db)
if __name__ == '__main__':
    getUpdateObjectPkl()