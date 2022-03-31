
features_to_index = {
    "RefAllel":5, 
    "VarAllel": 6, 
    "TotalPosCov":7, 
    "positioncoverage" : 8,
    "refForwardcoverage" : 9,
    "refReversecoverage" : 10,
    "varsCountOnForward" : 11, 
    "VarsCountOnReverse" : 12,
    "genotype" : 13,
    "frequency" : 14,
    "strandbiasflag" : 15, "meanPosition" : 16,
    "pstd" : 17,
    "meanQuality" : 18,
    "qstd" : 19, 
    "pvalue" : 20,
    "ratio" : 21,
    "mapq" : 22,
    "qratio" : 23,
    "higreq" : 24, 
    "extrafreq" : 25,
    "shift3" : 26,
    "msi" : 27,
    "msint" : 28,
    "nm" : 29,
    "hicnt" : 30,
    "hicov" : 31,
    "leftSequence" : 32,
    "rightSequence" : 33,
    "region" : 34,
    "varType" : 35,
    "duprate" : 36
    }

All_features = [
    "RefAllel", 
    "VarAllel", 
    "TotalPosCov", 
    "positioncoverage",
    "refForwardcoverage",
    "refReversecoverage",
    "varsCountOnForward", 
    "VarsCountOnReverse",
    "genotype",
    "frequency",
    "strandbiasflag",
    "meanPosition",
    "pstd",
    "meanQuality",
    "qstd", 
    "pvalue",
    "ratio",
    "mapq",
    "qratio",
    "higreq", 
    "extrafreq",
    "shift3",
    "msi",
    "msint",
    "nm",
    "hicnt",
    "hicov",
    "leftSequence",
    "rightSequence",
    "region",
    "varType",
    "duprate"
    ]
fvc_selected_features = [
     "TotalPosCov", 
     "positioncoverage",
     "refForwardcoverage",
     "refReversecoverage",
     "varsCountOnForward", 
     "VarsCountOnReverse",
    "frequency",
     "meanPosition",
     "pstd",
     "meanQuality",
     "pvalue",
     "ratio",
    "mapq",
     "qratio",
     "higreq", 
    "shift3",
    "msi",
     "msint",
    "nm",
     "hicnt",
     "hicov",
    ]

#som_features with fisher
som_features = [
    "Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt",
    "Var1Depth", "Var1AltDepth", "Var1RefFwdReads", "Var1RefRevReads", "Var1AltFwdReads", "Var1AltRevReads", "Var1Genotype", "Var1AF",
    "Var1Bias", "Var1PMean", "Var1PStd", "Var1QMean", "Var1QStd", "Var1MQ", "Var1Sig_Noise", "Var1HiAF", "Var1ExtraAF", "Var1NM", "Var1Pvalue", "Var1Oddr",
    "Var2Depth", "Var2AltDepth", "Var2RefFwdReads", "Var2RefRevReads", "Var2AltFwdReads", "Var2AltRevReads", "Var2Genotype", "Var2AF",
    "Var2Bias", "Var2PMean", "Var2PStd", "Var2QMean", "Var2QStd", "Var2MQ", "Var2Sig_Noise", "Var2HiAF", "Var2ExtraAF", "Var2NM", "Var2Pvalue", "Var2Oddr",
    "shift3", "MSI", "MSI_NT", "5pFlankSeq", "3pFlankSeq", "Seg", "VarLabel", "VarType",
    "Duprate1", "SV_info1", "Duprate2", "SV_info2", "Pvalue", "Oddratio"
    ]
som_features_to_index = {
    "Sample":0, "Gene":1, "Chr":2, "Start":3, "End":4, "Ref":5, "Alt":6,
    "Var1Depth":7, "Var1AltDepth":8, "Var1RefFwdReads":9, "Var1RefRevReads":10, "Var1AltFwdReads":11, "Var1AltRevReads":12, "Var1Genotype":13, "Var1AF":14,
    "Var1Bias":15, "Var1PMean":16, "Var1PStd":17, "Var1QMean":18, "Var1QStd":19, "Var1MQ":20, "Var1Sig_Noise":21, "Var1HiAF":22, "Var1ExtraAF":23, "Var1NM":24, "Var1Pvalue":25, "Var1Oddr":26,
    "Var2Depth":27, "Var2AltDepth":28, "Var2RefFwdReads":29, "Var2RefRevReads":30, "Var2AltFwdReads":31, "Var2AltRevReads":32, "Var2Genotype":33, "Var2AF":34,
    "Var2Bias":35, "Var2PMean":36, "Var2PStd":37, "Var2QMean":38, "Var2QStd":39, "Var2MQ":40, "Var2Sig_Noise":41, "Var2HiAF":42, "Var2ExtraAF":43, "Var2NM":44, "Var2Pvalue":45, "Var2Oddr":46,
    "shift3":47, "MSI":48, "MSI_NT":49, "5pFlankSeq":50, "3pFlankSeq":51, "Seg":52, "VarLabel":53, "VarType":54,
    "Duprate1":55, "SV_info1":56, "Duprate2":57, "SV_info2":58, "Pvalue":59, "Oddratio":60
    }

som_selected_features = [
    "Var1Depth", "Var1AltDepth", "Var1RefFwdReads", "Var1RefRevReads", "Var1AltFwdReads", "Var1AltRevReads", "Var1AF",
    "Var1PMean", "Var1PStd", "Var1QMean", "Var1QStd", "Var1MQ", "Var1Sig_Noise", "Var1HiAF", "Var1ExtraAF", "Var1NM", "Var1Pvalue", "Var1Oddr",
    "Var2Depth", "Var2AltDepth", "Var2RefFwdReads", "Var2RefRevReads", "Var2AltFwdReads", "Var2AltRevReads", "Var2AF",
    "Var2PMean", "Var2PStd", "Var2QMean", "Var2QStd", "Var2MQ", "Var2Sig_Noise", "Var2HiAF", "Var2ExtraAF", "Var2NM", "Var2Pvalue", "Var2Oddr",
    "shift3", "MSI", "MSI_NT",
    "Pvalue", "Oddratio"
    ]

som_rf_input_features = [
    "VarLabel", "VarType", "RefLength", "AltLength",
    "Var1DepthFreq", "Var1RefFwdFreq", "Var1AltFwdFreq",  "Var1AF",
    "Var1PMean", "Var1PStd", "Var1QMean", "Var1QStd", "Var1MQ", "Var1Sig_Noise", "Var1HiAF", "Var1ExtraAF", "Var1NM", "Var1Pvalue", "Var1Oddr",
    "Var2DepthFreq", "Var2RefFwdFreq", "Var2AltFwdFreq",  "Var2AF",
    "Var2PMean", "Var2PStd", "Var2QMean", "Var2QStd", "Var2MQ", "Var2Sig_Noise", "Var2HiAF", "Var2ExtraAF", "Var2NM", "Var2Pvalue", "Var2Oddr",
    "shift3", "MSI", "MSI_NT",
    "Pvalue", "Oddratio"
    ]

#som_rf_input_features = som_selected_features
