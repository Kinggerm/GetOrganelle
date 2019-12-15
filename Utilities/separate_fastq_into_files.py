#!/usr/bin/env python
import os
import sys
try:
    directories = set([os.path.join(sys.argv[1], x) for x in os.listdir(sys.argv[1]) if os.path.isdir(os.path.join(sys.argv[1], x))])
    for i in [x for x in os.listdir(sys.argv[1]) if x.endswith('fq.gz')]:
        this_dir = i.split('_')[0]
        if os.path.join(sys.argv[1], this_dir) not in directories:
            os.mkdir(os.path.join(sys.argv[1], this_dir))
            directories.add(os.path.join(sys.argv[1], this_dir))
        print(i+" > "+str(os.path.join(this_dir, i)))
        os.rename(os.path.join(sys.argv[1], i), os.path.join(sys.argv[1], this_dir, i))
except Exception as e:
    print(str(e))
    print("\nUsage: python separate.py dir")