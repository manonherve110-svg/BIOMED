import json
import requests

url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/property/CanonicalSMILES/JSON"


print("connection API ....")

reponse = requests.get(url)
data = reponse.json()

with open("data/pubchem_aspirin.json" , "w") as f:
    json.dump(data, f)

print("fin telechargement json")