
def get_bed_sorted(bed):
  regions = dict()
  with open(bed, 'r') as f:
    for line in f:
      chr, start, end = line.split('\t')[:3]
      if str(chr) in regions:
        regions[str(chr)].append([int(start), int(end)])
      else:
        regions[str(chr)] = list()
        regions[str(chr)].append([int(start), int(end)])
  for k, v in regions.items():
    regions[k] = sorted(v)
  return regions

def is_in_highconf(chr, start, end, regions):
  regs = regions[chr]
  #print(len(regs))
  #print('processin: ', chr + ":" + str(start) + ":" + str(end))
  for rstart, rend in regs:
    if end < rstart:
      break
    if start > rstart and start < rend:
      #print('in region:', start, rstart, rend)
      return True
  return False

def filter_by_region(data, highconf_regions):
  data['in_highconf'] = data.apply(lambda x: is_in_highconf(x['Chr'], x['Start'], x['End'], highconf_regions), axis = 1)
  return data[data['in_highconf'] == True]
  
def get_truth(truth_file):
  truth_vars = set()
  with open(truth_file, 'r') as f:
    for var in f:
      if var[0] == '#': continue
      items = var.split('\t')
      chrom, pos, id, ref, alt, _, filter = items[:7]         
      if (filter.find('PASS') == -1) and (filter.strip() != '.') :
        continue
      if len(chrom) < 6:
        site = chrom + ":" + pos + ":" + ref.upper() + ":" + alt.upper()
        truth_vars.add(site)
    return truth_vars

def get_label(truth_vars, key):
  if key in truth_vars:
    return 1
  else :
    return 0
  return 0
