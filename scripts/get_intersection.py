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
  #--case: insertion
  if len(ref) < len(alt):
      ref = "."
      alt = alt[1:]
  elif len(ref) > len(alt):
      start = str(int(start) + 1)
      alt = "."
      ref = ref[1:]
  #return(chrom, start, ref, alt, vaf)                
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

def main():
  vardict_vcf_file = "/home/haoz/data/tmpvcf.vcf"
  mutect_vcf_file = "/home/haoz/data/tmpvcf.vcf"
  strelka_vcf_file = "/home/haoz/data/tmpvcf.vcf"
  mc2v = {'mutect2':set(), 'vardict':set(), 'strelka2':set()} #map: caller to vairants set
  #map: caller to file name
  mc2f = {'mutect2':mutect_vcf_file, 'vardict':vardict_vcf_file, 'strelka2':strelka_vcf_file} 
  callers = mc2v.keys()
  for c in callers:
    mc2v[c] = get_variant_set(mc2f[c])

  for i in callers:
    line = ""
    for j in callers:
      if i == j: 
        line += "0\t" 
        continue
      #print("length of caller {} and {} set: {} {}".format(i, j, len(mc2v[i]), len(mc2v[j])))
      line += str(get_intersection_num(mc2v[i], mc2v[j]))
      line += '\t'
    line += '\n'
    print(line)

if __name__ == "__main__":
  main()
