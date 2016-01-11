#!/usr/bin/python
class Particle:
    
    IDUP    = '0'
    ISTUP   = '0'
    MOTHUP1 = '0'
    MOTHUP2 = '0'
    ICOLUP1 = '0' 
    ICOLUP2 = '0'
    PX      = '0'
    PY      = '0'
    PZ      = '0'
    M       = '0'
    VTIMUP  = '0'
    SPINUP  = '0'
    parent  =   0
    daughter1 = 0
    daughter2 = 0
    def __init__(self, line):
        self.IDUP    = line.split()[0]
        self.ISTUP   = line.split()[1]
        self.MOTHUP1 = line.split()[2]
        self.MOTHUP2 = line.split()[3]
        self.ICOLUP1 = line.split()[4]
        self.ICOLUP2 = line.split()[5]
        self.PX      = line.split()[6]
        self.PY      = line.split()[7]
        self.PX      = line.split()[8]
        self.E       = line.split()[9]
        self.M       = line.split()[10]
        self.VTIMUP  = line.split()[11]
        self.SPINUP  = line.split()[12]
        self.parent_ind  = 0
        self.daughter1_ind = 0
        self.daughter2_ind = 0
    def string(self):
        line = self.IDUP + '\t'
        line += self.ISTUP + '\t'
        line += '1\t'
        line += '2\t'
        line += str(self.parent_ind) + '\t'
        line += str(self.daughter1_ind) + '\t'
        line += str(self.daughter2_ind) + '\t'
        line += self.ICOLUP1 + '\t'
        line += self.ICOLUP2 + '\t'
        line += self.PX + ' '
        line += self.PY + ' '
        line += self.PZ + ' '
        line += self.E + ' '
        line += self.M + ' '
        line += self.VTIMUP + ' '
        line += self.SPINUP + ' '
        return line

    def string_orig(self):
        line = self.IDUP + '\t'
        line += self.ISTUP + '\t'
        line += str(int(self.MOTHUP1) - 1) + '\t'
        line += str(int(self.MOTHUP2) - 1) + '\t'
        line += self.ICOLUP1 + '\t'
        line += self.ICOLUP2 + '\t'
        line += str(self.parent_ind) + '\t'
        line += str(self.daughter1_ind) + '\t'
        line += str(self.daughter2_ind) + '\t'
        
        line += self.PX + '\t'
        line += self.PY + '\t'
        line += self.PZ + '\t'
        #line += self.E + ' '
        #line += self.M + ' '
        #line += self.VTIMUP + ' '
        #line += self.SPINUP + ' '
        return line


    def AddDaughter(self, daughter_ind):
        if self.daughter1_ind:
            self.daughter2_ind = daughter_ind
        else:
            self.daughter1_ind = daughter_ind
    
