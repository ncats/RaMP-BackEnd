from src.rampConfig.RampConfig import RampConfig
from src.util.EntityBuilder import EntityBuilder

resourceConfFile = "../config/external_resource_config.txt"
resourceConf = RampConfig()
resourceConf.loadConfig(resourceConfFile)

builder = EntityBuilder(resourceConf)
builder.download_pubchem_molecular_information()
