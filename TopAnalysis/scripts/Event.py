from Particle import *

class Event:
    Header        = ""
    Particles     = []

    def __init__(self):
        self.Header = ""
        self.Particles = []

    def ls(self):
        print self.Header + '\n'
        for ind in range(0, len(self.Particles)):
            line = self.Particles[ind].string_for_ls()
            print ind, "\t", line 

    def write(self, FILE):
        FILE.write('<event>\n')
        self.PrepareHeaderFor_write()
        FILE.write(self.Header + '\n')
        for ind in range(0, len(self.Particles)):
            line = self.Particles[ind].string_for_write()
            if line != "":
                FILE.write(line + '\n')
        
    def ConnectToParents(self, daughter_ind):
        daughter = self.Particles[daughter_ind]
        if (daughter.ISTUP == '-1'):
            self.Particles[daughter_ind].parent = 0
            return
        parent1 = self.Particles[int(daughter.MOTHUP1) -1]
        parent2 = self.Particles[int(daughter.MOTHUP2) -1]
        if daughter.MOTHUP1 != daughter.MOTHUP2:
            if parent1.ICOLUP1 == daughter.ICOLUP1 or parent1.ICOLUP2 == daughter.ICOLUP2:
                self.Particles[int(daughter.MOTHUP1) -1].AddDaughter(daughter_ind)
                self.Particles[daughter_ind].parent_ind = int(daughter.MOTHUP1) -1
            else:
                self.Particles[int(daughter.MOTHUP2) -1].AddDaughter(daughter_ind)
                self.Particles[daughter_ind].parent_ind = int(daughter.MOTHUP2) -1
        else:
            self.Particles[int(daughter.MOTHUP1) -1].AddDaughter(daughter_ind)
            self.Particles[daughter_ind].parent_ind = int(daughter.MOTHUP1) -1


    def Flip(self):
        for ind in range(0, len(self.Particles)) :
            self.ConnectToParents(ind)
        lepton_count = 0
        
        for ind in range(0, len(self.Particles)) :
            IDUP = self.Particles[ind].IDUP
            if IDUP == '24' or IDUP == '-24':
                Wboson = self.Particles[ind]
                daughter1 = self.Particles[Wboson.daughter1_ind]
                daughter2 = self.Particles[Wboson.daughter2_ind]
                if daughter1.ICOLUP1 == '0' and daughter1.ICOLUP2 == '0':
                    lepton_count += 1
                    continue
                parent             = self.Particles[Wboson.parent_ind]
                lightquark_ind     = 0
                antilightquark_ind = 0
                bquark_ind         = 0
                if daughter1.ICOLUP1 == '0':
                    lightquark_ind     = Wboson.daughter2_ind
                    antilightquark_ind = Wboson.daughter1_ind
                else:
                    lightquark_ind     = Wboson.daughter1_ind
                    antilightquark_ind = Wboson.daughter2_ind
                if self.Particles[parent.daughter1_ind].IDUP == '24' or self.Particles[parent.daughter1_ind].IDUP == '-24':
                    bquark_ind = parent.daughter2_ind
                else:
                    bquark_ind = parent.daughter1_ind

                if IDUP == '24':
                    Wboson_colour = parent.ICOLUP1
                    Wboson_anticolour = self.Particles[antilightquark_ind].ICOLUP2

                    self.Particles[bquark_ind].ICOLUP1 = Wboson_anticolour
                    self.Particles[lightquark_ind].ICOLUP1 = Wboson_colour
                else:
                    Wboson_colour = self.Particles[lightquark_ind].ICOLUP1
                    Wboson_anticolour = parent.ICOLUP2

                    self.Particles[bquark_ind].ICOLUP2 = Wboson_colour
                    self.Particles[antilightquark_ind].ICOLUP2 = Wboson_anticolour
                
                         
        return lepton_count
    def PrepareHeaderFor_write(self):
        N_INOUT_particles = 0
        for ind in range(0, len(self.Particles)) :
            if self.Particles[ind].ISTUP == '1' or self.Particles[ind].ISTUP == '-1':
                N_INOUT_particles += 1
        Header = str(N_INOUT_particles) + ' '
        self.Header = self.Header.lstrip()
        self.Header = Header + self.Header[len(self.Header.split()[0]):]
        
        
