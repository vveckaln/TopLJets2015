from Event import *
# This is a script for taking a ttbar LHE file and flipping the color strings between the W decay products and the b-quark from the top decay.
#
# contact: Benjamin Nachman (bnachman@cern.ch)
#

#flip=True

def WriteLHE():
    jlcount = 0
    jjcount = 0
    llcount = 0
    event_count = 0
    OUTPUTFILE = open('events_flipped_powheg_test.lhe', 'w')
    INPUTFILE = open('/afs/cern.ch/user/m/mseidel/public/TT_13TeV_powheg.lhe', 'r')
    event_count = 0
    event = 0
    for line in INPUTFILE:
        #print "line %s length %u"% (line, len(line.split()))
        start_line = False
        if '<event>' in line:
            event_count += 1
            event = Event()
            start_line = True
        if event:
            
            header_line = False
            particle_line = False
            if len(line.split()) == 6:
                event.Header = line
                header_line = True
            if len(line.split()) == 13:
                particle = Particle(line)
                event.Particles.append(particle)
                particle_line = True
            if not start_line and not header_line and not particle_line:
                print "   **** LISTING ORIGINAL **** "
                event.ls()
                leptons = event.Flip()
                
                print "   **** LISTING FLIPPED **** "
                event.ls()
                print " DONE DONE DONE DONE DONE DONE DONE DONE DONE "
 #               raw_input("Press Enter")
                if leptons == 0:
                    jjcount += 1
                elif leptons == 1:
                    jlcount += 1
                else:
                    llcount += 1
                event.write(OUTPUTFILE)
                event = 0
#            print start_line, header_line, particle_line
        if not event:
            OUTPUTFILE.write(line)

    OUTPUTFILE.write('</LesHouchesEvents>\n')
    print "semileptonic",   jlcount
    print "fully hadronic", jjcount
    print "dileptonic",     llcount




