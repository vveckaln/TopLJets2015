#!/usr/bin/env python2

from collections import OrderedDict

plots = OrderedDict()
plots["a"] = 1
for key in plots:
    print key

print plots["a"]
print "test ", plots["b"]
print "probe"
