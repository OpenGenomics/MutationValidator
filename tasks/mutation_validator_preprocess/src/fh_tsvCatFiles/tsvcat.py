#!/usr/bin/env python2.6
"""
tsvcat [files]

Concatenates TSV-with-header files, aligning columns with same name.  
Can rename columns and match columns across files with different names.

Modified by Lee Lichtenstein on June 20, 2012 to enable additional functionality and enable integration into Firehose
"""

import sys, itertools
import tsvutil
tsvutil.fix_stdio()
import csv
from OrderedSet import OrderedSet
#import codecs; sys.stdout = codecs.open('/dev/stdout','w',encoding='utf8',buffering=0)
import codecs

def flatten(iter):
    return list(itertools.chain(*iter))
  
def stable_uniq(x):
    s = set(); y = []
    for i in x:
        if i in s: continue
        s.add(i)
        y.append(i)
    return y
  
def tsv_reader(f):
    return csv.DictReader(f, dialect=None, delimiter="\t", quoting=csv.QUOTE_NONE)


def tsvcat(items, outFilePointer):
    alias_specs = [s for s in items if '=' in s]
    drop_specs = [s for s in items if s.startswith('-')]
    filenames = [s for s in items if s not in alias_specs and s not in drop_specs]
    #files = [open(f, 'r') for f in filenames]
    # codecs.open(filename, mode[, encoding

    files = [codecs.open(f, 'r', encoding="iso8859-1") for f in filenames]
    if not filenames:
        files = [sys.stdin]
    # LTL: Did not handle prepended comment lines.
    # file_cols = [f.readline()[:-1].split("\t") for f in files]
    # Unique comment lines are prepended onto the final result.
    file_cols = []
    commentLines = OrderedSet()
    for f in files:
        isHeaderFound = False
        line = f.readline()
        while not isHeaderFound:
            if not line.startswith('#'):
                isHeaderFound = True
                file_cols.append(line[:-1].split("\t"))
            else:
                commentLines.add(line.rstrip())
                line = f.readline()
    all_cols = stable_uniq(flatten(file_cols))
    
    aliases = {}
    for alias_spec in alias_specs:
        left, right = alias_spec.split('=')
        assert left != right
        assert left and right
        assert left in all_cols
        aliases[right] = left
        if right not in all_cols:
            all_cols[ all_cols.index(left) ] = right
        else:
            all_cols.remove(left)
    
    for drop_spec in drop_specs:
        col = drop_spec[1:]
        #print col
        #print all_cols
        assert col in all_cols
        all_cols.remove(col)
    
    if len(commentLines) > 0:
        outFilePointer.write("\n".join(commentLines)+"\n")

    outFilePointer.write("\t".join(all_cols)+"\n")
    
    for i, f in enumerate(files):
        cols = file_cols[i]
        for line in f:
            #line = unicode(line,'utf8')
            parts = line[:-1].split("\t")
            hash = {}
            for j in range(len(cols)):
                hash[cols[j]] = parts[j]
            out = []
            for col in all_cols:
                if col in hash:                
                    out.append(hash[col])
                elif col in aliases:
                    out.append(hash[aliases[col]])
                else:
                    out.append('')
            outFilePointer.write(u"\t".join(out) + '\n')


if __name__ == '__main__':
    items = sys.argv[1:]
    tsvcat(items,sys.stdout)

