import time
import sqlalchemy as sqla
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy import Column, Integer, String, ForeignKey, Float, Boolean,ForeignKey\
,ForeignKeyConstraint

'''
This script create a RaMP database with all tables based on the schema

'''
class RaMP_schema():
    # define all parameters for the database.
    password = 'Ehe131224'
    username = 'root'
    host = 'localhost'
    dbname = 'mathelabramp'
    engine = sqla.create_engine('mysql+pymysql://{}:{}@{}/{}'\
                                         .format(username,password,host,dbname))
        
    if not database_exists(engine.url):
        create_database(engine.url)
        
    
    Session = sessionmaker(bind = engine)
    session = Session()
    Base = declarative_base()
    Base.metadata.create_all(engine)


    class Analyte(Base):
        __tablename__ = 'analyte'
        rampId = Column(String(30),primary_key = True)
        type = Column(String(30))
        sources = relationship('Source',back_populates = 'analyte')
        
        def __repr__(self):
            return 'Analyte {}: Type {}'.format(self.rampId,self.type)
    
    class Analytehasontology(Base):
        __tablename__ = 'analytehasontology'
        rampCompoundId = Column(String(30),primary_key = True)
        rampOntologyIdLocation = Column(String(30),primary_key = True)
        def __repr__(self):
            return 'Analyte {}: Ontology: {}'.format(self.rampCompoundId,self.rampOntologyIdLocation)
    
    class Analytehaspathway(Base):
        __tablename__ = 'analytehaspathway'
        rampId = Column(String(30),primary_key = True)
        pathwayRampId = Column(String(30),primary_key = True)
        pathwaySource = Column(String(30))
        def __repr__(self):
            return 'Analyte {}: Pathway {}: PathwaySource {}'.format(self.rampId,
                                                                     self.pathwayRampId,
                                                                     self.pathwaySource)
    class Analytesynonym(Base):
        __tablename__ = 'analytesynonym'
        Synonym = Column(String(500),primary_key = True)
        rampId = Column(String(30),ForeignKey('analyte.rampId'))
        geneOrCompound = Column(String(30))
        source = Column(String(30))
        def __repr__(self):
            return 'Synonym {}: rampId {}: geneOrCompound {}: source {}'\
                .format(self.Synonym,self.rampId,self.geneOrCompound,self.source)
    
    class Catalyzed(Base):
        __tablename__ = 'catalyzed'
        rampCompoundId = Column(String(30),primary_key = True)
        rampGeneId = Column(String(30),primary_key = True)
        
        def __repr__(self):
            return 'Compound {}: Gene {}'\
                .format(self.rampCompoundId,self.rampGeneId)
    
    class Ontology(Base):
        __tablename__ = 'ontology'
        rampOntologyIdLocation = Column(String(30),primary_key = True)
        commonName = Column(String(30))
        biofluidORcellular = Column(String(30))
        def __repr__(self):
            return 'Ontology {}: CommonName {}: Type {}'\
                .format(self.rampOntologyIdLocation,
                        self.commonName,
                        self.biofluidORcellular)
    class Pathway(Base):
        __tablename__ = 'pathway'
        pathwayRampId = Column(String(30))
        sourceId =Column(String(30),primary_key = True)
        type = Column(String(30))
        pathwayCategory = Column(String(30))
        pathwayName = Column(String(30))
        
        def __repr__(self):
            return 'Pathway {}: Source Id {}: Type {}: Category {}: Name {}'\
                .format(self.pathwayRampId,self.sourceId,
                        self.type,self.pathwayCategory,self.pathwayName)
    
    class Source(Base):
        __tablename__ = 'source'
        sourceId = Column(String(500),primary_key = True)
        rampId = Column(String(30),ForeignKey('analyte.rampId'))
        IDtype = Column(String(30))
        geneOrCompound = Column(String(30))
        commonName = Column(String(30))
        analyte = relationship('Analyte',back_populates = 'sources')
        def __repr__(self):
            return '|SourceId {}: RampId {}: IDtype {}: geneOrCompound {}: commonName {}|'\
                .format(self.sourceId,self.rampId,self.IDtype,self.geneOrCompound,
                        self.commonName)
            
            
        


    
        