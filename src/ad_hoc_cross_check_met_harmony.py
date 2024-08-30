from src.rampConfig.RampConfig import RampConfig
from src.util.EntityBuilder import EntityBuilder

resourceConfFile = "../config/external_resource_config.txt"
resourceConf = RampConfig()
resourceConf.loadConfig(resourceConfFile)

builder = EntityBuilder(resourceConf)
builder.crossCheckMetaboliteHarmony(buildMetAndCompoundProps = True, criteria = "MW", tolerance = 0.1, pctOrAbs = 'pct')
