#!/usr/bin/env python3

def get_variants(fname)
  truth_indels = dict()
  filters = set()
  with open(truth_file, 'r') as f:
    for var in f:
      if var[0] == "#": continue
      items = var.split('\t')
      chrom, pos, id, ref, alt, _, filter = items[:7]         
      #if len(chrom) < 6 and filter == "PASS" and (len(ref) > 1 or len(alt) > 1) :
      #if filter.find("PASS") == -1 or filter.find("HighConf") == -1: continue
      filters.add(filter)
      if len(chrom) < 6 and (len(ref) > 1 or len(alt) > 1) :
        for alt_i in alt.split(","):
          site = chrom + ":" + pos + ":" + ref + ":" + alt_i
          truth_indels[site] = []

def get_prec_by_af(af_low, af_high):
  pass

if __name__ == "__main__":
  truth_file = ''
  truth_indels = get_variation(truth_indels)
  infile = ''
  rabv = get_variation(infile)
