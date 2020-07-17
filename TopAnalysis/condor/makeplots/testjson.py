#!/usr/bin/env python2
import os
from collections import OrderedDict

PROJECT=os.getenv("PROJECT")

import json
# c = 2, 1
# print c[1]
# # c[1] = 4

# a = []
# a += [[2, 0]]
# a += [[1, 4]]
# a += [[5, 6]]
# for i, j in a:
#     print i, j
# print a[1][1]
# print a

# a[1][1] = 9
# print a[1][1]
systdiction  = []
systjsonfilename = PROJECT + '/data/era2016/samples.json'
# #jsonfilename = 'newtestjson.json'
systjsonFile = open(systjsonfilename, 'r')

systdiction += json.load(systjsonFile, encoding='utf-8', object_pairs_hook=OrderedDict).items()

print systdiction
print "-------------"
print systdiction[0]
print systdiction["MC13TeV_SingleTbar_tW"]
#print diction
# , encoding='utf-8', object_pairs_hook=OrderedDict).items()
# systjsonFile.close()

# samplejsonfilename = PROJECT + '/data/era2016/samples.json'
# samplejsonFile = open(samplejsonfilename, 'r')
# samplediction = json.load(samplejsonFile)
# samplejsonFile.close()
# samplediction = {"sample0" : [5, "b"], "sample1" : [6, "c"]}
# templist = []
# sample = []
# for syst in systdiction:
#     a = systdiction[syst]
#     print a
#     templist += [a, "syst", False] 
#     print templist
#     sample += syst
# for slist, sample, booltry in templist:
#     print slist
#     print "\n"

