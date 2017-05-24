'''
Created on 16 Jun 2013

@author: stewart
'''

import argparse
import sys
import os
import re
from tsvcat import tsvcat


if not (sys.version_info[0] == 2  and sys.version_info[1] in [7]):
    raise "Must use Python 2.7.x"


def parseOptions():
    description = '''
    Filter a tsv file based on fields conditions from command line arguments
    '''

    epilog= '''

        Writes subset of input lines satisfying field condition(s)

        -i <id> -m <input tsv> -f <list of fields> -a <list of binary operations> -v <list of values> -c <list of connectors> -o <output area>

            '''
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('-i','--id', metavar='id', type=str, help ='id.')
    parser.add_argument('-t','--tsvFiles', metavar='tsvFiles', type=str, help ='input tsv files.')
    parser.add_argument('-c','--criteria', metavar='criteria', type=str, help ='selection criteria of form (<field> [compare] [value]) [op] ... .')
    parser.add_argument('-o','--output', metavar='output', type=str, help ='output area.',default='.')
    parser.add_argument('-k','--keys', metavar='keys', type=str, help ='unique key fields.')
    parser.add_argument('-x','--extension', metavar='extension', type=str, help ='file name extension.',default='tsv')
    #parser.add_argument('-d','--delimiter', metavar='delimiter', type=str, help='delimiter',default='\t')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parseOptions()
    id = args.id
    tsvFiles = args.tsvFiles.split(',')
    criteria = args.criteria
    keys = args.keys.split(',')
    extension = args.extension
    #delimiter = args.delimiter
    #if "\\" in delimiter:
    #    delimiter = ""
    output = args.output

    print "command line:"
    print " ".join(sys.argv) + "\n"

    # parse criteria for fields in <brackets>
    fields = ''
    if not criteria is None:
        fields=re.findall('\[\s?(.+?)\s?\]', criteria)

    NT=len(tsvFiles)
    NF=len(fields)
    kf= [-1] * NF

    if output is None:
        output = "."

    if not os.path.exists(output):
        os.makedirs(output)

    print "input tsvFiles:\n\t" + "\n\t".join(tsvFiles) + "\n"
    print "input keys:\n\t" + "&".join(keys) + "\n"


    output_filename = output + "/" + id + "." + extension
    print "output file:\t" + output_filename + "\n"
    output_filenames= ['.'] * NT

    nf = -1
    for tsvFile in tsvFiles:
        nf = nf + 1
        # loop over tsv file
        output_filenames[nf]=output_filename + '.' + str(nf) + '.tsv'
        outputFileFP = file(output_filenames[nf], 'w')
        fh = open(tsvFile,'rt')
        nin = 0
        npass = 0
        nreject = 0
        n = 0
        nhead = 0

        for lines in fh:
            n = n+1
            line = lines.strip()
            if (n == 1) & (line[0] == '#') & (not criteria is None):
                outputFileFP.write(line + ": Filtered: " + criteria + "\n")
                continue

            if (line[0] == '#'):
                outputFileFP.write(line+ "\n")
                continue

            x = line.split()

            if (nhead==0):
                i=0
                for i in range(NF):
                    f=fields[i]
                    if x.count(f)>0:
                        kf[i]=x.index(f)
                    else:
                        print "missing field "+f
                        continue



                if any([k1<0 for k1 in kf]):
                    raise NameError('tsv file missing required fields ')

                outputFileFP.write(line+"\n")
                nhead=1
                continue

            nin=nin+1
            KEEP = True
            if not criteria is None:
                C = criteria.replace('[','"').replace(']','"')

                for i in range(NF):
                    x1=x[kf[i]]
                    C=C.replace(fields[i],x[kf[i]])

                KEEP = eval(C)

            if KEEP:
                outputFileFP.write(line + "\n")
                npass=npass+1

            else:
                nreject=nreject+1


        fh.close()
        outputFileFP.close()

        print "\ninput  file:\t" + tsvFile
        print "input  #calls:\t" + str(nin)
        print "pass   #calls:\t" + str(npass)
        print "reject #calls:\t" + str(nreject) + "\n"

    if not (keys is None):
        output_filename1 = output + "/" + id +  '.pre-uniquified.' + extension
    else:
        output_filename1 = output_filename

    outputFileFP = file(output_filename1, 'w')
    tsvcat(output_filenames, outputFileFP)
    outputFileFP.close()
    print "output  file:\t" + output_filename

    outputFileFP = file(output_filename, 'w')


    if not (keys is None):
        NK=len(keys)
        kf= [-1] * NK
        KEYLIST=[]

        fh = open(output_filename1,'rt')
        nin = 0
        npass = 0
        nreject = 0
        n = 0
        nhead = 0

        for lines in fh:
            n = n+1
            line = lines.strip()
            if (line[0] == '#'):
                outputFileFP.write(line+ "\n")
                continue

            x = line.split()

            if (nhead==0):
                i=0
                for i in range(NK):
                    f=keys[i]
                    if x.count(f)>0:
                        kf[i]=x.index(f)
                    else:
                        print "missing field "+f
                        continue



                if any([k1<0 for k1 in kf]):
                    raise NameError('tsv file missing required key fields for uniquifying')

                nhead=1
                outputFileFP.write(line+"\n")
                continue

            nin=nin+1

            xk=[x[k1] for k1 in kf]
            key = '&'.join(xk)
            if not (key in KEYLIST):
                outputFileFP.write(line+"\n")
                npass=npass+1;
            else:
                nreject=nreject+1


        fh.close()
        outputFileFP.close()

        print "\nraw  file:\t" + output_filename1
        print "\nuniquified  file:\t" + output_filename

        print "input  #calls:\t" + str(nin)
        print "pass   #calls:\t" + str(npass)
        print "reject #calls:\t" + str(nreject) + "\n"



