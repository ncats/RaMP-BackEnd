from updateSQL import RampUpdater
import pickle as pk
if __name__ == '__main__':
    pkf = open('../misc/output/wikipathwayRdfPk.pkl','rb')
    wp = pk.load(pkf)
    rp = RampUpdater(wp)
    rp.checkNewEntry()