#!/usr/bin/env python3

#############################################################################
# This module is designed for get the consensus. such as the vardict and 
# mutect2 
#############################################################################

def get_info_vardict(var):
  #Chrom_     Start_  End     Ref_    Alt_    Gene    Type    VAF_    CHGVS   PHGVS
  #chrom start end ref alt gene vaf
  info = var[7].split(";")
  vaf = info[5]
  if vaf[:3] == "AF=":
      vaf = float(vaf[3:])
  else:
      print("AF position error! are you sure???")
      exit(0)
  chrom = ""
  if(var[0][:3] == "chr"):
      chrom = var[0][3:]
  else:
      chrom = var[0]
  #TODO start and end??
  var_type = info[2]
  if var_type[:4] == "TYPE": var_type = var_type[5:]
  start = 0
  ref = ""
  alt = ""
  if var_type == "Insertion":
      start = var[1]
      ref = "."
      alt = var[4][1:]
  elif var_type == "Deletion":
      start = str(int(var[1]) + 1)
      ref = var[3][1:]
      alt = "."
  else:
      start = var[1]
      ref = var[3]
      alt = var[4]
  return tuple((chrom, start, ref,
                alt, vaf))

NOINDEL=False
NOSNV=False
def get_vcf_info(var):
  #chrom start end ref alt gene vaf
  var = var.split('\t')
  chrom = ""
  if(var[0][:3] == "chr"):
      chrom = var[0][3:]
  else:
      chrom = var[0]
  info = var[7].split(":")
  vaf = "0.1"
  start = var[1]
  ref = var[3]
  alt = var[4]
  #snv:
  if NOSNV and len(ref) == 1 and len(alt) == 1:
    return ""

  if NOINDEL and (len(ref) > 1 or len(alt) > 1):
    return ""
  
  #--case: insertion
  # if len(ref) < len(alt):
  #     ref = "."
  #     alt = alt[1:]
  # elif len(ref) > len(alt):
  #     start = str(int(start) + 1)
  #     alt = "."
  #     ref = ref[1:]
  #return(chrom, start, ref, alt, vaf)                
  #print("debug:", chrom+":"+start+":"+ref+":"+alt)
  return chrom+":"+start+":"+ref+":"+alt

def get_variant_set(filename):
  variant_set = set()
  with open(filename, 'r') as f:
    for line in f:
      if line[0] == '#': continue
      variant_set.add(get_vcf_info(line))
  return variant_set

def get_intersection_num(s1, s2):
  counter = 0
  for i in s1:
    if i in s2:
      counter += 1
  return counter

def padding(s):
  return s + ''.join((2-len(s)//8)*'\t')

def main():
  #vardict_vcf_file = "/home/user_home/haoz/workspace/RabbitVar/benchresults/RabVar_T75N05.vcf"
  #mutect_vcf_file = "/home/user_home/haoz/bcbio_bench_result/FDSyn_impure/T75N05_work/final_mutect2/SPP_GT_3-1_1/SPP_GT_3-1_1-mutect2.vcf"
  #strelka_vcf_file = "/home/user_home/haoz/bcbio_bench_result/FDSyn_impure/T75N05_work/final_strelka2/SPP_GT_3-1_1/SPP_GT_3-1_1-strelka2.vcf"
  vardict_vcf_file = "/home/user_home/haoz/workspace/RabbitVar/benchresults/RabVar_FD2.vcf"
  mutect_vcf_file = "/home/user_home/haoz/workspace/Experiments/variant_calling/FD_T2-mutect2.vcf"
  strelka_vcf_file = "/home/user_home/haoz/workspace/Experiments/variant_calling/FD_T2-strelka2.vcf"
  ground_truth_file = "/home/user_home/haoz/data/FD/Truth/FDALL_TRUTH.vcf"
  mc2v = {'mutect2':set(), 'vardict':set(), 'strelka2':set(), 'gt':set()} #map: caller to vairants set
  #map: caller to file name
  mc2f = {'mutect2':mutect_vcf_file, 'vardict':vardict_vcf_file, 'strelka2':strelka_vcf_file, 'gt':ground_truth_file} 
  callers = mc2v.keys()
  for c in callers:
    mc2v[c] = get_variant_set(mc2f[c])

  for i in callers:
    line = padding(i)
    base = len(mc2v[i])
    for j in callers:
      if i == j: 
        line += "{}\t".format(len(mc2v[i]) / base) 
        continue
      #print("length of caller {} and {} set: {} {}".format(i, j, len(mc2v[i]), len(mc2v[j])))
      line += str(get_intersection_num(mc2v[i], mc2v[j]) / base)
      line += '\t'
    line += '\n'
    print(line)

if __name__ == "__main__":
  main()
