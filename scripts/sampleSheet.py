#! /usr/bin/env python3
from sys import argv
import re

final = []
with open(argv[1], 'rt') as IN, open(argv[2], 'wt') as OUT:
    text = IN.read()
    lines = text.split('\n')
    index = [i for i, word in enumerate(lines) if word.startswith('Sample_Name')][0] 
    final = lines[0:index+1]
    for L in lines[index:]:
        F = L.split(',')
        m = re.match("^" + argv[3] + "$", F[0])
        if m != None:
            final.append(L)

    OUT.write('\n'.join(final))
