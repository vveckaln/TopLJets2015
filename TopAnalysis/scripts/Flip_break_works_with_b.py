# This is a script for taking a ttbar LHE file and flipping the color strings between the W decay products and the b-quark from the top decay.
#
# contact: Benjamin Nachman (bnachman@cern.ch)
#

#flip=True

def unsplit(A):
    out=''
    for a in A:
        out=out+' '
        out=out+a
        pass
    return out

def islepton(A):
    if abs(A)>10 and abs(A)<17:
        return True
    else:
        return False
    pass

def WriteLHE(flip,inDir):


    A=0
    B=0
    C=0

    oecount=0
    elcount=0
    llcount=0
    jjcount=0

    import os
    count=0
    fileList=[]
    if os.path.isfile(inDir) : fileList.append(inDir)
    else : fileList=[f for f in os.listdir(inDir) if os.path.isfile(os.path.join(inDir, f))]

    os.system('mkdir -p outfiles/')
    os.system('rm outfiles/*')

    for a in fileList :
        count+=1

        print 'Flipping events in ',a
        F2 = open('events_flipped_powheg_%d.lhe'%count, 'w')
        F = open(a,'r')

        #print oecount,begin,target
        begin=0
        target=999999999999999
        if (oecount-begin > target):
            break
        
        if (count==1):
            #First, write the header.
            for f in F:
                if '<event>' in f:
                    break
                F2.write(f)
                pass
            pass

        event=[]
        ecount=0
        for f in F:
            if (elcount > target):
                break
            if '<event>' in f:
                oecount+=1
                #print "     ",oecount,begin,target
                if (oecount - begin > target):
                    break
                if (ecount> 0 and oecount > begin):
                    topnum=0
                    #F2.write('<event>\n')
                    #First, remove the tops and W's.  To do this we need to first reduce the particle counter at the top
                    #Note: we could keep the leptonic W to get tau polarizations correct.  For later. 
                    counter=1
                    W1=-1
                    W2=-1
                    for ei in event:
                        if (len(ei.split())==13):
                            if ei.split()[0]=='-24' or ei.split()[0]=='-6' or ei.split()[0]=='24' or ei.split()[0]=='6':
                                topnum+=1
                                pass
                            if (ei.split()[0])=='24':
                                W1=counter
                                pass
                            if (ei.split()[0])=='-24':
                                W2=counter
                                pass
                            counter+=1
                            pass
                        pass
                    #We only want semi-leptonic ttbar
                    lepcounter=0
                    saveW=-1
                    hadW=-1
                    for ei in event:
                        if (len(ei.split())<10):
                            continue
                        if abs(int(ei.split()[0]))==11 or abs(int(ei.split()[0]))==13 or abs(int(ei.split()[0]))==15:
                            lepcounter+=1
                            if int(ei.split()[2])==W1:
                                saveW=W1
                                hadW=W2
                            else:
                                saveW=W2
                                hadW=W1
                                pass
                            pass
                        pass
                    #print lepcounter,W1,W2,saveW,hadW
                    mycounter=1
                    fi=[]
                    li=[]
                    bcolor=-1
                    wcolor=-1
                    if (lepcounter==0):
                        #F2.write('<event>\n')
                        jjcount+=1
                        for ei in event:
                            if 'event' in ei:
                                continue
                            #F2.write(ei)
                            pass
                        pass
                    if (lepcounter==2):
                        F2.write('<event>\n')
                        llcount+=1
                        for ei in event:
                            if 'event' in ei:
                                continue
                            F2.write(ei)
                            pass
                        pass
                    if (lepcounter > 2):
                        print "error, more than 2 leptons"
                        pass
                    if (lepcounter==1):
                        F2.write('<event>\n')
                        pass
                    for ei in event:
                        myskip=False
                        if (lepcounter!=1):
                            break
                        #print ei
                        if (len(ei.split())==6):
                            ei=str(int(ei.split()[0])-topnum)+' '+unsplit(ei.split()[1:len(ei.split())])+'\n'
                            F2.write(ei)
                            pass
                        if (len(ei.split())<10):
                            continue
                        if ei.split()[0]=='5' and ei.split()[1]!='-1' and ei.split()[2]!='1':
                            bcolor=int(ei.split()[4])
                            fi+=[ei]
                            myskip=True
                            pass
                        if ei.split()[2]==str(hadW):
                            wcolor=max(int(ei.split()[4]),int(ei.split()[5]))
                            fi+=[ei]
                            myskip=True
                            pass
                        if ei.split()[1]!='-1' and (abs(int(ei.split()[0]))<=16 or abs(int(ei.split()[0])==24)):
                            ei=unsplit([ei.split()[0],ei.split()[1],str(1),str(2)]+ei.split()[4:len(ei.split())])+'\n'
                            pass
                        if mycounter==saveW:
                            #F2.write(ei)
                            fi+=[ei]
                        if ei.split()[0]!='-6' and ei.split()[0]!='-24' and ei.split()[0]!='6' and ei.split()[0]!='24' :
                            #print "write"
                            #print ei
                            #if (islepton(int(ei.split()[0]))):
                            #    li+=[ei]
                            #else:
                            if (myskip):
                                pass
                            else:
                                F2.write(ei)
                                fi+=[ei]
                                pass
                            pass
                        mycounter+=1
                        pass
                    mycounter=1
                    saveW=-1
                    fi=fi+li
                    #print 'bcolor',bcolor,'wcolor',wcolor
                    #Now, we go back through, find the new position of the leptonic W and set the lepton mother to this W. 
                    for ei in fi:
                        if len(ei.split())<10:
                            continue
                        if abs(int(ei.split()[0]))==24:
                            saveW=mycounter
                            pass
                        if islepton(int(ei.split()[0])):
                            ei=unsplit([ei.split()[0],ei.split()[1],str(saveW),str(saveW)]+ei.split()[4:len(ei.split())])+'\n'
                            #F2.write(ei)
                            pass
                        mycounter+=1
                        pass
                    #Now for the flipping
                    for ei in fi:
                        if len(ei.split())<10:
                            continue
                        if ei.split()[0]=='5' and ei.split()[1]!='-1' and ei.split()[2]!='1':
                            if (flip):
                                ei=unsplit([ei.split()[0],ei.split()[1],str(1),str(2),str(wcolor),str(0)]+ei.split()[6:len(ei.split())])+'\n'
                            else:
                                ei=unsplit([ei.split()[0],ei.split()[1],str(1),str(2),str(bcolor),str(0)]+ei.split()[6:len(ei.split())])+'\n'
                                pass
                            F2.write(ei)
                            elcount+=1
                            pass
                        if ei.split()[2]==str(hadW):
                            if (int(ei.split()[0]))<0:
                                ei=unsplit([ei.split()[0],ei.split()[1],str(1),str(2),str(0),str(wcolor)]+ei.split()[6:len(ei.split())])+'\n'
                                pass
                            else:
                                if (flip):
                                    ei=unsplit([ei.split()[0],ei.split()[1],str(1),str(2),str(bcolor),str(0)]+ei.split()[6:len(ei.split())])+'\n'
                                else:
                                    ei=unsplit([ei.split()[0],ei.split()[1],str(1),str(2),str(wcolor),str(0)]+ei.split()[6:len(ei.split())])+'\n'
                                    pass
                                pass
                            F2.write(ei)
                            pass
                        pass
                    
                    if (lepcounter>0):
                        F2.write('</event>\n') #+str(elcount+llcount)+'\n')
                        pass
                    pass
                    
                #if ecount>2:
                #    break
                event=[]
                ecount+=1
                pass
            event+=[f]

        F2.write('</LesHouchesEvents>\n')

    A+=elcount
    B+=jjcount
    C+=llcount
    print "semileptonic",elcount
    print "fully hadronic",jjcount
    print "dileptonic",llcount




