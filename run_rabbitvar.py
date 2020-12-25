import joblib
import numpy as np
import pandas as pd
import time 
import argparse
import subprocess
import pybedtools

def get_socks():
    for line in str(subprocess.check_output('lscpu')).split('\\n'):
        contex = line.split()
        if contex[0].strip() == 'NUMA':
            return contex[-1]
def prepare_cmd(bed_file_name, out_file_name):
    cmd = '''{BIN} \
    -i /home/old_home/haoz/workspace/data/NA12878/{BED} \
    -G /home/old_home/haoz/workspace/data/hg38/hg38.fa \
    -f 0.01 \
    -N FD_T_2|FD_N_2 \
    -b {FDDir}/FD_T2.sorted.bam|{FDDir}/FD_N2.sorted.bam  \
    -c 1 -S 2  -E 3 -g 4   \
    --th 20 --fisher \
    --auto_resize \
    --out {OUT} \
    '''.format(BIN="/home/old_home/haoz/workspace/FastVC/build/FastVC",
               BED=bed_file_name,
               FDDir='/home/ssd',
               OUT=out_file_name)
    cmd = cmd.split()
    return cmd


if __name__ == "__main__":
    socks = get_socks()
    print(socks)
    #re-distribure bed file
    beds = pybedtools.example_bedtool('/home/old_home/haoz/workspace/data/NA12878/mini1.bed')
    sorted(beds)
    nregions = len(beds)
    #beds[0:nregions//2].saveas('mini1.bed')
    #beds[nregions//2: ].saveas('mini2.bed')
    ps = list()
    for i in range(int(socks)):
        cmd = prepare_cmd('mini{}.bed'.format(i+1),
                          'tmp{}.txt'.format(i+1))
        ps.append(subprocess.Popen(cmd, stderr=subprocess.DEVNULL))

    for p in ps:
        p.wait()

    print("now, all process run over!")
