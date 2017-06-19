
import csv
import sys

input_maf = sys.argv[1]
output_maf = sys.argv[2]
maf_dialect = sys.argv[3]


ifid = open(input_maf)
ofid = open(output_maf,'w')

line = ifid.readline()
while line.startswith('#'):
    ofid.write(line)
    line = ifid.readline()

if maf_dialect == 'TCGA':
    line.replace('Start_position','Start_Position')
    line.replace('End_position','End_Position')
elif maf_dialect == 'Broad':
    line.replace('Start_Position','Start_position')
    line.replace('End_Position','End_position')
else:
    raise Exception('maf_dialect must be TCGA or Broad')

while line:
    ofid.write(line)
    line = ifid.readline()

