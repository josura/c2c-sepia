import numpy as np
import pandas as pd



metalinks = pd.read_csv('/home/josura/Downloads/name_map.csv')
metalinks.to_csv('metalinks.csv')

metabolites = pd.read_csv('/home/josura/Projects/ccc/fluxes/scFEA/data/cName_complete_mouse_c70_m168.csv')

metalinks_kegg = metalinks[metalinks["KEGG"].notnull()]["KEGG"]

metabolites_kegg = list(metabolites.iloc[0,:])

metabolites_fluxes = pd.read_csv('/home/josura/Downloads/name_map(1).csv')

#dropping na  for the HMDB column

metabolites_fluxes = metabolites_fluxes[metabolites_fluxes["HMDB"].notnull()]

metalinks_hmdb = metalinks[metalinks["HMDB"].notnull()]["HMDB"]
metabolites_hmdb = list(metabolites_fluxes["HMDB"])

#set intersection between first row of metabolites and metalinks["KEGG"]
common = list(set(metabolites_kegg) & set(metalinks_kegg))
common_hmdb = list(set(metabolites_hmdb) & set(metalinks_hmdb))

common_pd = metalinks[metalinks["KEGG"].isin(common)]
common_pd.to_csv('common_kegg.csv')