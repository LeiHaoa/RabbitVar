#	str.emplace_back(sample);
#	str.emplace_back(region.gene);
#	str.emplace_back(region.chr);
#       str.emplace_back(std::to_string(variant->startPosition));
#       str.emplace_back(std::to_string(variant->endPosition));
#       str.emplace_back(variant->refallele);
#       str.emplace_back(variant->varallele);
#       str.emplace_back(std::to_string(variant->totalPosCoverage));
#       str.emplace_back(std::to_string(variant->positionCoverage));
#       str.emplace_back(std::to_string(variant->refForwardCoverage));
#       str.emplace_back(std::to_string(variant->refReverseCoverage));
#       str.emplace_back(std::to_string(variant->varsCountOnForward));
#       str.emplace_back(std::to_string(variant->varsCountOnReverse));
#       str.emplace_back(variant->genotype == "" ? "0" : variant->genotype);
#       str.emplace_back(std::to_string(variant->frequency));
#       str.emplace_back(variant->strandBiasFlag);
#       str.emplace_back(std::to_string(variant->meanPosition));
#       str.emplace_back(variant->isAtLeastAt2Positions ? std::to_string(1) :std::to_string(0));//ptsd
#       str.emplace_back(std::to_string(variant->meanQuality)); //qual
#       str.emplace_back(variant->hasAtLeast2DiffQualities ? std::to_string(1) : std::to_string(0)); //qstd
#       str.emplace_back(std::to_string(variant->meanMappingQuality)); //mapq
#       str.emplace_back(std::to_string(variant->highQualityToLowQualityRatio)); //qratio 
#       str.emplace_back(std::to_string(variant->highQualityReadsFrequency)); //higreq
#       str.emplace_back(std::to_string(variant->extraFrequency)); //extrafreq
#       str.emplace_back(std::to_string(variant->shift3)); //shift3
#       str.emplace_back(std::to_string(variant->msi)); //msi
#       str.emplace_back(std::to_string(variant->msint)); //msint
#       str.emplace_back(variant->numberOfMismatches > 0 ? std::to_string(variant->numberOfMismatches) : std::to_string(0)); // nm
#       str.emplace_back(std::to_string(variant->hicnt)); //hicnt
#       str.emplace_back(std::to_string(variant->hicov)); //hicov
#       str.emplace_back(variant->leftseq.empty() ? "0": variant->leftseq); //leftSequence
#       str.emplace_back(variant->rightseq.empty() ? "0": variant->rightseq); //rightSequence
#       str.emplace_back(region.chr + ":" + std::to_string(region.start) + "-" + std::to_string(region.end)); //region
#       str.emplace_back(variant->vartype); //varType
#       str.emplace_back(std::to_string(variant->duprate)); //duprate
#       str.emplace_back(sv == "" ? "0": sv); //sv
#       str.emplace_back(std::to_string(variant->crispr));
#compare c and java result data 

#sample_name chr1    chr1    3829690 3829690 T   C   43  43  0   0   17  26  C/C 1.0000  0;2 31.0    1
#(meanquality)38.6    1   60.0
#(qratio)86.000  1.0000  0
#(sift3)0   2.000   2   1.1
#(hicht)43  43  ACGCCGGGGCTTCCCATACA    GTGGCTCCTTCAAATGACAT    chr1:3829690-3918526
#SNV 0   0
def is_same_var(jri, cri):
    
    if(    jri[3] == cri[3] #start position
       and jri[4] == cri[4] #end position   
       and jri[5] == cri[5] #refallele      
       and jri[6] == cri[6] #varallele      
       and int(jri[7]) == int(cri[7]) #totalposcoverage
       and int(jri[8]) == int(cri[8]) #positioncoverage
       and int(jri[9]) == int(cri[9]) #refForwardcoverage
       and int(jri[10]) == int(cri[10]) #refReversecoverage
       and int(jri[11]) == int(cri[11]) #varsCountOnForward
       and int(jri[12]) == int(cri[12]) #VarsCountOnReverse
       and jri[13] == cri[13] #genotype
       and abs(float(jri[14]) - float(cri[14])) <= 0.001 #frequency
       and jri[15] == cri[15] #strandbiasflag
       and abs(float(jri[16]) - float(cri[16])) <= 0.1 #meanPosition
       and jri[17] == cri[17] #pstd
       and abs(float(jri[18]) - float(cri[18])) <= 0.1 #meanQuality 
       and jri[19] == cri[19] #qstd
       and abs(float(jri[20]) - float(cri[20])) <= 0.1 #mapq
       and abs(float(jri[21]) - float(cri[21])) <= 0.01 #qratio
       and abs(float(jri[22]) - float(cri[22])) <= 0.001 #higreq
       and abs(float(jri[23]) - float(cri[23])) <= 0.001 #extrafreq
       and int(jri[24]) == int(jri[24]) #shift3
       and abs(float(jri[25]) - float(cri[25])) <= 0.01 #msi
       and int(jri[26]) == int(cri[26]) #msint
       and abs(float(jri[27]) - float(cri[27])) <= 0.1 #nm
       and int(jri[28]) == int(cri[28]) #hicnt
       and int(jri[29]) == int(cri[29]) #hicov
       and jri[30] == cri[30]
       and jri[31] == cri[31]
       and jri[32] == cri[32]
       and jri[33] == cri[33]
    ):
        return True
 
    else:
        return False

def is_same_var_stat(jri, cri):
    same_count = 0
    if jri[3] == cri[3]: #start position
        same_count += 1 
    else:
        return False
    if jri[4] == cri[4]: #end position
        same_count += 1 
    else:
        return False
    if jri[5] == cri[5] : #refallele      
        same_count += 1 
    if jri[6] == cri[6] : #varallele      
        same_count += 1 
    if int(jri[7]) == int(cri[7]) : #totalposcoverage
        same_count += 1 
    if int(jri[8]) == int(cri[8]) : #positioncoverage
        same_count += 1 
    if int(jri[9]) == int(cri[9]) : #refForwardcoverage
        same_count += 1 
    if int(jri[10]) == int(cri[10]) : #refReversecoverage
        same_count += 1 
    if int(jri[11]) == int(cri[11]) : #varsCountOnForward
        same_count += 1 
    if int(jri[12]) == int(cri[12]) : #VarsCountOnReverse
        same_count += 1 
    if jri[13] == cri[13] : #genotype
        same_count += 1 
    if abs(float(jri[14]) - float(cri[14])) <= 0.001 : #frequency
        same_count += 1 
    if jri[15] == cri[15]: #strandbiasflag
        same_count += 1 
    if abs(float(jri[16]) - float(cri[16])) <= 0.1 : #meanPosition
        same_count += 1 
    if jri[17] == cri[17] :#pstd
        same_count += 1 
    if abs(float(jri[18]) - float(cri[18])) <= 0.1 : #meanQuality 
        same_count += 1 
       #and jri[18] == cri[18]
    if jri[19] == cri[19]: #qstd
        same_count += 1
    if abs(float(jri[20]) - float(cri[20])) <= 0.1: #mapq
        same_count += 1
    if abs(float(jri[21]) - float(cri[21])) <= 0.01: #qratio
        same_count += 1
    if abs(float(jri[22]) - float(cri[22])) <= 0.001: #higreq
        same_count += 1
    if abs(float(jri[23]) - float(cri[23])) <= 0.001: #extrafreq
        same_count += 1
    if int(jri[24]) == int(jri[24]): #shift3
        same_count += 1
    if abs(float(jri[25]) - float(cri[25])) <= 0.01: #msi
        same_count += 1
    if int(jri[26]) == int(cri[26]): #msint
        same_count += 1
    if abs(float(jri[27]) - float(cri[27])) <= 0.1: #nm
        same_count += 1
    if int(jri[28]) == int(cri[28]): #hicnt
        same_count += 1
    if int(jri[29]) == int(cri[29]): #hicov
        same_count += 1
    if jri[30] == cri[30]:
        same_count += 1
    if jri[31] == cri[31]:
        same_count += 1
    if jri[32] == cri[32]:
        same_count += 1
    if jri[33] == cri[33]:
        same_count += 1

    if same_count >= 28:
        #if same_count != 31:
        #    print(same_count);
        #    print("|".join(jri) ,"\n", "|".join(cri))
        return True
    else:
        return False
    
def compare(cpp_result, java_resutl):
    cr = list()
    jr = list()
    #cr_len = list()
    #jr_len = list()
    with open(cpp_result, 'r') as f:
        for var in f:
            items = var.split('\t')
            #cr_len.append(len(items))
            if(len(items) > 35):
                cr.append(items)
        
    with open(java_resutl, 'r') as f:
        for var in f:
            items = var.split('\t')
            #jr_len.append(len(items))
            if(len(items) > 35):
                jr.append(items)

    cpp_len = len(cr)
    java_len = len(jr)
    result = 0
    for jr_i in jr:
        for cr_i in cr:
            #print("comparing: ", "|".join(jr_i), "|".join(cr_i))
            if is_same_var_stat(jr_i, cr_i):
            #if is_same_var(jr_i, cr_i):
                result += 1
                break
    #print("len info: \n", cr_len, jr_len)
    print("accurcy: ", cpp_len, java_len, result, result/float(cpp_len))


if __name__ == "__main__":
    cpp_result = "./tmp.vcf"
    java_result = "../VarDictJava/tmp.vcf"
    compare(cpp_result, java_result)
    
    #500 region result:
    #('accurcy: ', 131793, 131800, 131661, 0.9989984293551251)
    #real    24m11.519s
