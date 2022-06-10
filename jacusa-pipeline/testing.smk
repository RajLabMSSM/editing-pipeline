#RNA Editing Pipline
#Winston Cuddleston, Raj Lab
#2022

import pandas as pd
import os

#references
refDir = config["refDir"]
print(refDir)
editingRefDir = config["editingRefDir"]
print(editingRefDir)
humandbDir = config["humandbDir"]
print(humandbDir)

#data
projectDir = config["projectDir"]
print(projectDir)
metadata =  pd.read_csv(config["metadata"], sep = "\t")
samples = metadata['sample']
print(samples)
metadata_dict = metadata.set_index('sample').T.to_dict()
print(metadata_dict)
