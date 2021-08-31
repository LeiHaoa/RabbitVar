#!/usr/bin/env python3

# This script is judge if a variant in a bed file regions
#     example: if an variant is an exon varaint

import multiprocessing as mp

def get_bed(bed_file):
  bed = dict()
  with open(bed_file, 'r') as bf:
    for line in bf:
      line = line.split("\t")
      chrm, start, end = line[0:3]
      if chrm in bed:
        bed[chrm].append([int(start), int(end)])
      else:
        bed[chrm] = list()
        bed[chrm].append([int(start), int(end)])

    for k, v in bed.items():
      bed[k] = sorted(v) 
  return bed

def judge_in(t, bed):
  chrm, start = t
  start = int(start)
  region_list = bed[chrm]
  find = False
  for region in region_list:
    s, e = region
    if start < s:
      break
    else: # start >= s
      if start <= e:
        find = True
        break
  return find

def subfunc(variants, bed, q):
  total = len(variants)
  exon = 0
  for item in variants:
    #chrm, start = item
    if judge_in(item, bed):
      exon += 1
  q.put([total, exon])
  #if judge_in(chrm, int(start), end, bed):
  #  exon += 1
def main(bed_file, vcf_file, proc_num):
  bed = get_bed(bed_file)
  g_exon = 0
  g_total = 0
  variants = []
  with open(vcf_file, 'r') as vf:
    for line in vf:
      if line[0] == "#":
        continue
      line = line.split("\t")
      chrm, start = line[0:2]
      variants.append([chrm, start])
  q = mp.Queue(proc_num)
  procs = []
  stride = int(len(variants) / proc_num)
  for p in range(proc_num):
    proc = mp.Process(target = subfunc, 
        args = (variants[p*stride:min((p + 1)*stride, len(variants))], bed,  q, )
        )
    procs.append(proc)
    proc.start()
  
  for p in procs:
    p.join()

  for p in range(proc_num):
    total, exon = q.get()
    g_total += total
    g_exon += exon

  print("totall {} variants, {} extrons, and {} introns".format(g_total, g_exon, g_total - g_exon))

if __name__ == "__main__":
  #bed_file = "/home/old_home/haoz/workspace/data/Twist_Exome_Target_hg38.bed"
  #bed_file = "/home/data/haoz/FD/trainingSet_0118/genome.bed"
  bed_file = "/home/old_home/haoz/workspace/VCTools/bcbio_nextgen/genomes/Hsapiens/hg38/coverage/capture_regions/Exome-NGv3.bed"
  #vcf_file = "/home/data/haoz/FD/Truth/FDtruth_Data_2.indel.vcf"
  #vcf_file = "/home/old_home/haoz/workspace/FastVC/detection_result/FD_DATASET/FD_DATA_2_VardictFormated.vcf"
  vcf_file = "/home/old_home/haoz/workspace/FastVC/RandomForest/RandomForest_Filtered.vcf"
  proc_num = 40
  main(bed_file, vcf_file, proc_num)
