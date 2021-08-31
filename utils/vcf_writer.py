################################################################ 
# this model contained functions to write record into VCF file # 
################################################################ 
import math
import re
from .Global import label_to_varLabel, label_to_types

def write_header(fout):
  header = """\
##fileformat=VCFv4.3
##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name (with whitespace translated to underscores)">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type: SNV Insertion Deletion Complex">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=END,Number=1,Type=Integer,Description="Chr End Position">
##INFO=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=BIAS,Number=1,Type=String,Description="Strand Bias Info">
##INFO=<ID=REFBIAS,Number=1,Type=String,Description="Reference depth by strand">
##INFO=<ID=VARBIAS,Number=1,Type=String,Description="Variant depth by strand">
##INFO=<ID=PMEAN,Number=1,Type=Float,Description="Mean position in reads">
##INFO=<ID=PSTD,Number=1,Type=Float,Description="Position STD in reads">
##INFO=<ID=QUAL,Number=1,Type=Float,Description="Mean quality score in reads">
##INFO=<ID=QSTD,Number=1,Type=Float,Description="Quality score STD in reads">
##INFO=<ID=SBF,Number=1,Type=Float,Description="Strand Bias Fisher p-value">
##INFO=<ID=ODDRATIO,Number=1,Type=Float,Description="Strand Bias Odds ratio">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mean Mapping Quality">
##INFO=<ID=SN,Number=1,Type=Float,Description="Signal to noise">
##INFO=<ID=HIAF,Number=1,Type=Float,Description="Allele frequency using only high quality bases">
##INFO=<ID=ADJAF,Number=1,Type=Float,Description="Adjusted AF for indels due to local realignment">
##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description="No. of bases to be shifted to 3 prime for deletions due to alternative alignment">
##INFO=<ID=MSI,Number=1,Type=Float,Description="MicroSatellite. > 1 indicates MSI">
##INFO=<ID=MSILEN,Number=1,Type=Float,Description="MicroSatellite unit length in bp">
##INFO=<ID=NM,Number=1,Type=Float,Description="Mean mismatches in reads">
##INFO=<ID=LSEQ,Number=1,Type=String,Description="5' flanking seq">
##INFO=<ID=RSEQ,Number=1,Type=String,Description="3' flanking seq">
##INFO=<ID=GDAMP,Number=1,Type=Integer,Description="No. of amplicons supporting variant">
##INFO=<ID=TLAMP,Number=1,Type=Integer,Description="Total of amplicons covering variant">
##INFO=<ID=NCAMP,Number=1,Type=Integer,Description="No. of amplicons don't work">
##INFO=<ID=AMPFLAG,Number=1,Type=Integer,Description="Top variant in amplicons don't match">
##INFO=<ID=HICNT,Number=1,Type=Integer,Description="High quality variant reads">
##INFO=<ID=HICOV,Number=1,Type=Integer,Description="High quality total reads">
##INFO=<ID=SPLITREAD,Number=1,Type=Integer,Description="No. of split reads supporting SV">
##INFO=<ID=SPANPAIR,Number=1,Type=Integer,Description="No. of pairs supporting SV">
##INFO=<ID=DUPRATE,Number=1,Type=Float,Description="Duplication rate in fraction">
##FILTER=<ID=Bias,Description="Strand Bias">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference forward, reverse reads">
##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="Variant forward, reverse reads">
"""
  fout.write(header)

def format_record(record):
  try:
    [sample, gene, chrt, start, end, ref, alt, dp1, vd1, rfwd1, rrev1, vfwd1, \
     vrev1, gt1, af1, bias1, pmean1, pstd1, qual1, qstd1, mapq1, sn1, hiaf1, \
     adjaf1, nm1, sbf1, oddratio1, dp2, vd2, rfwd2, rrev2, vfwd2, vrev2, \
     gt2, af2, bias2, pmean2, pstd2, qual2, qstd2, mapq2, sn2, hiaf2, \
     adjaf2, nm2, sbf2, oddratio2, shift3, msi, msilen, lseq, rseq, seg, \
     status, vtype, duprate1, sv1, duprate2, sv2, pvalue, oddratio]  = record[:61]
  except ValueError:
    print("invalide record: \n", record, "\n record length ---> ", len(record))
    exit(-1)
  
  rd1 = rfwd1 + rrev1
  rd2 = rfwd2 + rrev2
  if vtype == "" : 
    vtype = "REF"
  GTFREQ = 0.2
  FREQ = 0.02
  gt  = "1/1" if (1-af1 < GTFREQ) else ("1/0" if af1 >= 0.5 else ("0/1" if af1 >= FREQ else "0/0"))
  gtm = "1/1" if (1-af2 < GTFREQ) else ("1/0" if af2 >= 0.5 else ("0/1" if af2 >= FREQ else "0/0"))
  bias1 = bias1.replace(';',',')
  bias2 = bias2.replace(';',',')
  if bias1 == "0": bias1 = "0,0" 
  if bias2 == "0": bias2 = "0,0" 
  mapq1 = f'{mapq1:.0}'
  mapq2 = f'{mapq2:.0}'
  qual = int(math.log(vd1)/math.log(2) * qual1) if vd1 > vd2  else int(math.log(vd2)/math.log(2) * qual2)

  pinfo1 = "\t".join([chrt, str(start), ".", ref, alt, str(qual)])
  filters = "PASS"
  sample_nowhitespace = re.sub(r'\s', '_', sample)

  pinfo2_1 = "STATUS={};SAMPLE={};TYPE={};DP={};VD={};AF={:.6f};SHIFT3={};MSI={};MSILEN={};SSF={};SOR={};LSEQ={};RSEQ={}".format(label_to_varLabel[status], sample_nowhitespace, label_to_types[vtype], dp1, vd1, af1, shift3, msi, msilen, pvalue, oddratio, lseq, rseq)
  pinfo2_2 = "{}:{}:{}:{}:{}:{}:{:.6f}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{:.6f}:{:.6f}:{}".format(gt, dp1, vd1, str(vfwd1) + "," + str(vrev1), str(rfwd1)+","+str(rrev1), str(rd1)+","+str(vd1), af1, bias1, pmean1, pstd1, qual1, qstd1, sbf1, oddratio1, mapq1, sn1, hiaf1, adjaf1, nm1) 
  pinfo2_3 = "{}:{}:{}:{}:{}:{}:{:.6f}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{:.6f}:{:.6f}:{}".format(gtm, dp2, vd2, str(vfwd2) + "," + str(vrev2), str(rfwd2)+","+str(rrev2), str(rd2)+","+str(vd2), af2, bias2, pmean2, pstd2, qual2, qstd2, sbf2, oddratio2, mapq2, sn2, hiaf2, adjaf2, nm2)
  pinfo2 = "\t".join([pinfo2_1, "GT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM", pinfo2_2, pinfo2_3])
  return "\t".join([pinfo1, filters, pinfo2])
