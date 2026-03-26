#This script only changes binning of approach = 1;
#From terminal run: python changeApproach1.py config-map.json <desired_binning_binPhi,(int)> <desired_binning_binZ,(int)> 
#Returns: config-map-binned.json

import sys
import json

if len(sys.argv) < 4:
    conFile = "config-map.json"
    bins = {"binPhi":sys.argv[1],"binZ":sys.argv[2]}
else:
    conFile = sys.argv[1]

with open(conFile, "r") as con_file:
    config = json.load(con_file)
    bins = {"binPhi":sys.argv[2],"binZ":sys.argv[3]}

surfaces = config["Surfaces"] #The surfaces dictionary containing keys "1", "2"...

for surface in surfaces:
        for object in surfaces[surface]: #An object in the list 'surface' contains a dictionary of properties
            if "approach" in object and object["approach"] == 1:
                 object["value"]["material"]["mapMaterial"] = True
                 #object["value"]["material"]["mappingType"] = "PreMapping" # IA
                 for bin_property in object["value"]["material"]["binUtility"]["binningdata"]:
                    if bin_property["value"] == "binPhi":
                        bin_property["bins"] = int(bins["binPhi"])
                    else: bin_property["bins"] = int(bins["binZ"])

for surface in surfaces:
        for object in surfaces[surface]: #An object in the list 'surface' contains a dictionary of properties
            if object["volume"]== 14: #10: # beam pipe!
                 object["value"]["material"]["mapMaterial"] = True
                 #object["value"]["material"]["mappingType"] = "PreMapping" # IA
                 for bin_property in object["value"]["material"]["binUtility"]["binningdata"]:
                    if bin_property["value"] == "binPhi":
                        bin_property["bins"] = int(bins["binPhi"])
                    else: bin_property["bins"] = int(bins["binZ"])


config["Surfaces"] = surfaces
con_file = json.dumps(config, indent=4)
with open('config-map-binned.json', 'w') as file:
    file.write(con_file)
