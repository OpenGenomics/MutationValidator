'''
Created on Jun 20, 2012

@author: lichtens
'''
import subprocess
import argparse
import sys
from tsvcat import tsvcat

def call(command, isPrintCmd=True):
    ''' returns returncode, output 
    '''
    try:
        if isPrintCmd:
            print command
        return 0, subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError as cpe:
        return cpe.returncode,cpe.output

def parseOptions():
    description = '''Given a tsv text file with headers on the first line, prune lines based on values seen in the specified headers.  The first column is ignored and the second is passed into tsvcat (which is used to concatenate the tsv files)'''
    epilog= '''
    This script passes a list of files through to tsvcat.py.  
    
    This script is necessary to get a large number of files processed in Firehose
    '''
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tsvFileListIn', metavar='tsvFileListIn',  type=str, help ='tsv file name in the Firehose format.  This is the file with a list of files to concatenate in the second column.  No headers nor comments in the tsv are accepted.')
    parser.add_argument('outputFilename', metavar='outputFilename',  type=str, help ='tsv file name for output.')
    
    args = parser.parse_args()
    
    return args 

if __name__ == '__main__':
    args = parseOptions()
    tsvFileListIn = args.tsvFileListIn
    outputFilename = args.outputFilename
    
    fpIn = file(tsvFileListIn, 'r')
    # fpOut = file(outputFilename, 'w')
    import codecs
    fpOut = codecs.open(outputFilename, 'w', encoding="iso-8859-1")
    
    fileList = []
    line = fpIn.readline()
    while line != '':
        if line.strip() != '':
            fileList.append(line.strip().split('\t')[1])
        line = fpIn.readline()
        pass
    
    if len(fileList) != 0:
        tsvcat(fileList, fpOut)
    pass
