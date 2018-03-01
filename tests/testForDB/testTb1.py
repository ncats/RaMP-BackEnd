import time
import csv
import unittest
import random
import time
from schema import RaMP_schema
from sqlalchemy import update, exc
class TestDatabase(unittest.TestCase):
    '''
    This test case prepend all IDtype to source ID
    so they are in format IDTYPE:SOURCEID
    '''
    def testMain(self):
        sche = RaMP_schema()
        sess = sche.session
        Source = sche.Source
        prepended_ids = sess.query(Source).filter(~Source.sourceId.like('%:%'))
        len(prepended_ids.all()
        for each in prepended_ids.all():
            print('Updating {} ...'.format(each))
            if ' ' not in each.sourceId:
                if each.IDtype == 'enzymeNomenclature':
                    each.sourceId = 'EN' + ':'+each.sourceId
                else:
                    each.sourceId = each.IDtype + ':'+each.sourceId.replace(' ','')
            else:
                sess.delete(each)
                
        sess.commit()
        
        
        

if __name__ == '__main__':
    unittest.main()