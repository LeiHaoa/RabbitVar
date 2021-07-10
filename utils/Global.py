

SOM_INDEL_FEATURES =  43 #46+3
SOM_SNV_FEATURES = 41 + 7 #41 + 7

#GERM_SNV_FEATYRES
#GERM_INDEL_FEATYRES

label_to_varLabel = {
   0 :  "Germline"      ,
   1 :  "StrongLOH"     ,
   2 :  "LikelyLOH"     ,
   3 :  "StrongSomatic" ,
   4 :  "LikelySomatic" ,
   5 :  "AFDiff"        ,
   6 :  "SampleSpecific",
}

label_to_types = {
        0: "SNV",
        1: "Deletion", 
        2: "Insertion", 
        3: "Complex", 
}
