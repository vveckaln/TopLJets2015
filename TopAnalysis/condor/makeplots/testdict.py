#!/usr/bin/env python2
import ROOT

if False:
    diff = 5.0
print diff
exit(0)
stra = "stra"
print "stra [" + stra + "]"

print "afdas",
print " asds"
exit(0)

vara = ROOT.TH1F("h", "h", 5, 0.0, 5.0)
print vara["key"]
exit(0)
d={}
print d==None, len(d)

string1 = "dsfds"
string2 = "string2"
print "printing something", string1, string2
exit(0)
diction = {
    "a0" : {
        "b0" : 1
    },
    "a1" : {
        "b1" : 1
    }
} 

diction3 = {
    "a30" : {
        "b30" : 1
    },
    "a31" : {
        "b31" : 1
    }
} 
diction2= {
    "key1" : {
        "key1_a2" : {
            "key1_b0" : 8,
            "key1_b1" : 4,

        "key1_b2" : 3
        }
    },
    "key2" :  {
        "key2_a2" : {
            "key2_b0" : 0,
            
        "key2_b2" : 3
        }
    }
}

diction4 = {"a" :0, "b":1, "c":2}
for key, numb in diction4:
    print key
    print numb
    print "===="
greatlist = [(diction, 0)]
greatlist = greatlist + [(diction3, 0)]
print greatlist
raw_input("penter")
for key in diction2:
    greatlist += [(diction2[key], 1)]
for dicti, val in greatlist:
    print dicti
    print val
    print "---"

#, (diction2, 1))

print greatlist

print diction

print diction2

#diction.update(diction2)
diction.update({"newkey": 10})
for key in diction.keys():
    if key == "a0":
        del diction[key]
    print key
# if diction.get("a1") == None: 
#     diction["a1"] = {}
# diction["a1"]["b0"] = 2
# print diction

# for s in (for s2 in diction[s]):
#     print s

# l = list(diction.keys())
# print l
print diction
exit(0)
