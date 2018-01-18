from update.RaMPDatabase import RaMPDatabase

class RaMPFixer(RaMPDatabase):
    '''
    This class simulates the function of the C# code in RaMP
    Get the original RaMP data from MySQL, then do the following things:
    1) For each table, find the missing data cell, and remove the entire row since the map
    is not successful. (Some columns are deliberated to left empty with 'NA'.)
    2) Remove special character such as ' ', and wrong data in the cell. Consider drop entire
    row based on conditions
    3) Remapping corrected RaMP ID relations based on corrected data.
    4) Create and import new RaMP data to the database.
    '''
    
    def __init__(self):
        super().__init__()
        
        
    
    def create_new_db(self):
        conn = self.connectToRaMP(dbname = None)
        
        with conn.cursor() as cur:
            query = "create database mathelabramp2;"
            cur.execute(query)
                    
            
        conn.commit()
            
        
        conn.close()
        
        print("hi")