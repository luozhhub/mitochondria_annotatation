"""
this file is request for mitochondria database:https://www.hmtvar.uniba.it/
"""
import requests
import re
import pandas as pd

def annotation():
    main_url = "https://www.hmtvar.uniba.it/api/main/"
    page_locus = "position"
    page_variant = "mutation"
    
    
    human_mitochondria_length = 16569
    pos_dict_list = []
    for i in range(1, human_mitochondria_length + 1):
        locus_response = requests.get(main_url + page_locus + "/" + str(i))
        if locus_response.status_code != 200:
            print("wrong in site: %s" % i)
        #print(locus_response.text)
        result = locus_response.json()
        
        for entry in result:
            ref = entry["ref_rCRS"]
            alt = entry["alt"]
            aa_change = entry["aa_change"]
            if aa_change is None：
                aa_1 = None
                aa_2 = None
                aa_pos = None
            else:
                aa = re.sub(r"\d", "", aa_change)
                aa_1 = aa[0]
                aa_2 = aa[1]
                aa_pos = re.sub(r"\D", "", aa_change)
            pos_dict_list.append({"position": i, "ref": ref, "alt": alt, "aa_old": aa_1, "aa_new": aa_2, "aa_pos": aa_pos})
            
    df = pd.DataFrame(pos_dict_list)
    print(df[0:5])
    df = df[["position", "ref", "alt", "aa_old", "aa_new", "aa_pos"]]
    df.to_csv("human_MT_mutation_annotation.csv", sep=",")
    
if __name__ == "__main__":
    annotation()
