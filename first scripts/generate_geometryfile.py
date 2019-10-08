# TODO:
#
#  x applyMovement routines accepting functions to create positions, e.g. sinus (in addition to linear/log)
#  => FOR NOW: implemented sin and sinrev as stepmode options, function seems not straightforward :/
#
#  x routine which plots the path in 3D
#  => IMPLEMENTED
#
#  x change loops that generate movement in applyMovement routines
#     -> make internal temporary list that holds all new positions, 
#        in the end just add this as a whole to the main lists
#        (no more appending and calculating further on previously appended item, simplifies 'addup' mode)
#        (easier to introduce input functions to determine movement)
#  => IMPLEMENTED
#  
#  x implement loop changes for log mode in campos
#  => IMPLEMENTED + log functionality now correct (log+logrev mode)
#
#  x implement circle routine
#  => IMPLEMENTED
#     x implement centering
#     x check lefthandedness
#
#  x implement helix routine
#  => IMPLEMENTED
#
#  x implement previous loop changes to lookat
#  => IMPLEMENTED
#  x implement addup mode in lookat (sky?)
#  => IMPLEMENTED
#
#  x plotPath also plots lookat?
#  x plotPath also plots sky?
#
#  x unequal axes in circle/helix
#  - implement radiusMovement method
#  - implement spiral, using radiusMovement and circle
#  - implement "pause" method
#  - implement campos/lookat loop mechanism to sky/fidx?
#  - 'sin' or 'function' mode: also give as input a range which should be mapped by the function to the movement
#    (e.g. for sin, give the x-range it should go through on it's path, so 0..Pi/2 or something)
#
#  - make helix/circle start from last known point, not starting circle always at a fixed point


import sys
import os
import math



class pathGenerator:
    """class to generate an animation path and write it to a geometry file to use in Splotch"""

    def __init__(self, filename = "test", campos = [10,0,0], lookat = [0,0,0], sky = [0,0,1], fidx = 0):

        self.campos = [campos]
        self.lookat = [lookat]
        self.sky = [sky]
        self.fidx = [fidx]

        self.firstCampos = campos
        self.firstLookat = lookat
        self.firstSky = sky
        self.firstFidx = fidx
        self.filename = filename

        self.lastAction = False

    def lastCampos(self):
        try:
            self.campos[-1][0]
            return list(self.campos[-1])
        except:
            self.campos = [self.firstCampos]
            return list(self.firstCampos)

    def lastLookat(self):
        try:
            self.lookat[-1][0]
            return list(self.lookat[-1])
        except:
            self.lookat = [self.firstLookat]
            return list(self.firstLookat)

    def lastSky(self):
        try:
            self.sky[-1][0]
            return list(self.sky[-1])
        except:
            self.sky = [self.firstSky]
            return list(self.firstSky)

    def lastFidx(self):
        try:
            self.fidx[-1]
            return self.fidx[-1]
        except:
            self.fidx = [self.firstFidx]
            return self.firstFidx

    def addSimpleLineSegment(self, startCampos, endCampos, nsteps, mode="linear", lookatFixed=True):
        """ simple straight line trajectory, fixed lookat and sky, fixed fidx """
        nBefore = len(self.campos)
        if startCampos == "continue":
            startCampos = self.lastCampos()

        self.addCamposMovement(startCampos, endCampos, nsteps, mode, applymode="expand")
        self.padLookat(nsteps)
        self.padSky(nsteps)
        self.padFidx(nsteps)

        if not lookatFixed: 
            endlookat = [lastlookat + (endcampos - startcampos) for (lastlookat,endcampos,startcampos) in zip(self.lastLookat(), endCampos, startCampos)]
            self.addLookatMovement(self.lastLookat(), endlookat, nsteps=nsteps, applymode="modify")

        self.lastAction = dict()
        self.lastAction["action"] = "addLineSegment"
        self.lastAction["nBefore"] = nBefore
        self.lastAction["nAfter"] = len(self.campos)
        self.lastAction["expandedLists"] = ["campos", "lookat", "sky", "fidx"]
        self.lastAction["modifiedLists"] = []

    def addCircle(self, radius, axis, rotation, direction, nsteps, center=[0,0,0], quiet=False):

        if not quiet:
            print "* adding circle to path [ r=",radius,"| ax=",axis,"| rot=",rotation,\
            "| center=",center,"]"

        nBefore = len(self.campos)

        try:
            len(radius)
        except:
            radius=[radius,radius]

        if axis=="x":
            axlist=[1,2]
        elif axis=="y":
            axlist=[2,0]
        elif axis=="z":
            axlist=[0,1]

        start = [0,0,0]
        start[axlist[0]] = radius[0]
        if self.lastCampos() != start:
            self.campos.append(start)
            nsteps = nsteps-1

        nquart = int(rotation*4)
        nquartSteps = int(nsteps/nquart)

        if direction=="righthanded":
            initsign=1
        elif direction=="lefthanded":
            initsign=-1
        sign = initsign

        sinDict = dict()
        sinDict[1] = "sinrev"
        sinDict[-1] = "sin"


        for q in range(1,nquart+1):
            end1 = self.lastCampos()
            end1[axlist[0]] += radius[0]*((-1)**(math.ceil(q*1./2)))
            end2=[0,0,0]
            end2[axlist[1]] += initsign*radius[1]*((-1)**(math.ceil(q*1./2+3./2)))
            self.addCamposMovement(start=self.lastCampos(), end=end1, nsteps=nquartSteps,
                                   stepmode=sinDict[initsign*sign], applymode="expand", quiet=True)
            self.addCamposMovement(start=[0,0,0], end=end2, nsteps=nquartSteps,
                                   stepmode=sinDict[initsign*(-sign)], applymode="addup", quiet=True)
            sign = -sign

        self.lastAction = dict()
        self.lastAction["action"] = "addCircle"
        self.lastAction["nBefore"] = nBefore
        self.lastAction["nAfter"] = len(self.campos)
        self.lastAction["modifiedLists"] = []
        self.lastAction["expandedLists"] = ["campos"]

        if center!=[0,0,0]:
            self.addCamposOffset(center, quiet=True)

        self.padLists()

    def addHelix(self, radius, axis, rotation, direction, nsteps, displacement, center=[0,0,0],
                 lookatFixed=True, stepmode="linear", quiet=False):

        if not quiet:
            print "* adding helix [ r=",radius,"| ax=",axis,"| rot=",rotation,\
            "| disp=",displacement,stepmode,"| center=",center,"]"

        nBefore = len(self.campos)

        if axis=="x":
            displacementAxis=0
        elif axis=="y":
            displacementAxis=1
        elif axis=="z":
            displacementAxis=2
        displacementStart=[0,0,0]
        displacementEnd=[0,0,0]
        displacementEnd[displacementAxis]=displacement

        self.addCircle(radius, axis, rotation, direction, nsteps, center=center, quiet=True)
        self.addCamposMovement(start=displacementStart, end=displacementEnd, nsteps=0,
                               stepmode=stepmode, applymode="addup", quiet=True)
        if not lookatFixed:
            self.addLookatMovement(start=displacementStart, end=displacementEnd, nsteps=0,
                                   stepmode=stepmode, applymode="modify", quiet=True)
            
        self.lastAction = dict()
        self.lastAction["action"] = "addHelix"
        self.lastAction["nBefore"] = nBefore
        self.lastAction["nAfter"] = len(self.campos)
        self.lastAction["modifiedLists"] = []
        self.lastAction["expandedLists"] = ["campos"]
        if not lookatFixed:
            self.lastAction["expandedLists"].append("lookat")

    def addCamposMovement(self, start, end, nsteps, stepmode="linear", applymode="expand", quiet=False):
        """ 
        basic routine to add movement in camera position
        - stepmode = 'linear', 'log', 'sin' or 'sinrev'
        - applymode = 'expand', 'modify' or 'addup'
          'expand' just adds new camera positions
          'modify' applies camera movement to previous action (if exists), possibly overwriting previously 
          generated camera positions (so input variable nsteps is not used)
          'addup' also works on the previous action (like 'modify'), but adds up the values to the existing 
          ones of the previous action, instead of deleting the previous ones. multiple addups are possible
        """

        if applymode!="expand" and self.lastAction!=False:
            index = self.lastAction["nBefore"]
            addupStorage = self.eraseList(index, self.campos)
            nsteps = self.lastAction["nAfter"] - self.lastAction["nBefore"]
            if applymode == "modify":
                if not quiet:
                    print "  * adding campos movement to previous action - deleting old positions"
                addupStorage = [[0,0,0] for i in range(0,nsteps)]
            elif applymode == "addup":
                if not quiet:
                    print "  * adding campos movement to previous action - adding up to old positions"
        else:
            addupStorage = [[0,0,0] for i in range(0,nsteps)]

        nBefore = len(self.campos)

        if self.lastCampos() != start and applymode != "addup":
            self.campos.append(start)
            nsteps = nsteps-1

        appendArray = []
        # ---------------------------------------------------------------------------------------------
        if stepmode=="linear":
            stepsizes = [(end[i]-start[i])*1./nsteps for i in range(0,3)]
            appendArray = [[start[i] + j*stepsizes[i] for i in range(0,3)] for j in range(1,nsteps+1)]
        # ---------------------------------------------------------------------------------------------
        elif stepmode=="log" or stepmode=="logrev":
            # some dirty tweaking here to have something well behaved (don't bother figuring out why and what :p)
            stepsizes = [math.pow(10,(math.log10(0.001)-math.log10(1.001))*1./nsteps)
                         if (start[i] != end[i]) else 1 for i in range(0,3)]
            if stepmode=="log":
                steplist = [[1-(stepsizes[i]**j - 0.001) if (start[i] != end[i]) else 1 for i in range(0,3)] 
                            for j in range(1,nsteps+1)]
            elif stepmode=="logrev":
                steplist = [[stepsizes[i]**j - 0.001 if (start[i] != end[i]) else 1 for i in range(0,3)] 
                            for j in range(0,nsteps)]
                steplist.reverse()
            steplist[-1] = [1.,1.,1.]
            appendArray = [[(end[i]-start[i])*steplist[j][i] for i in range(0,3)] 
                           for j in range(0,nsteps)]
        # ---------------------------------------------------------------------------------------------
        elif stepmode=="sin" or stepmode=="sinrev":
            amplitudes = [(end[i]-start[i])*1. for i in range(0,3)]
            stepsize = math.pi/(2*nsteps)
            xarray = [j*stepsize for j in range(1,nsteps+1)]
            if stepmode=="sin":
                yarray = [math.sin(x) for x in xarray]
            elif stepmode=="sinrev":
                yarray = [1-math.cos(x) for x in xarray]
            for j in range(0,len(yarray)):
                appendArray.append([start[i] + yarray[j]*amplitudes[i] for i in range(0,3)]) 
        # ---------------------------------------------------------------------------------------------         
        else:
            print "# ERROR : mode not recognized ["+mode+"]"
            exit(-1)
        # ---------------------------------------------------------------------------------------------

        if applymode=="addup":
            appendArray = [[addupStorage[j][i] + appendArray[j][i] for i in range(0,3)] for j in range(0,nsteps)]

        self.campos.extend(appendArray)

        self.lastAction = dict()
        self.lastAction["action"] = "addCamposMovement"
        self.lastAction["nBefore"] = nBefore
        self.lastAction["nAfter"] = len(self.campos)
        self.lastAction["expandedLists"] = []
        self.lastAction["modifiedLists"] = []
        if applymode == "expand":
            self.lastAction["expandedLists"] = ["campos"]
        else:
            self.lastAction["modifiedLists"] = ["campos"]

    def addLookatMovement(self, start, end, nsteps, stepmode="linear", applymode="expand", quiet=False):
        """ 
        basic routine to add movement in lookat position
        - stepmode = 'linear' or 'log'
        - applymode = 'expand' or 'modify'
          'expand' just adds new lookat positions
          'modify' applies lookat movement to previous action (if exists), possibly overwriting previously 
          generated lookat positions (so input variable nsteps is not used)
        """

        if applymode!="expand" and self.lastAction!=False:
            index = self.lastAction["nBefore"]
            addupStorage = self.eraseList(index, self.lookat)
            nsteps = self.lastAction["nAfter"] - self.lastAction["nBefore"]
            if applymode == "modify":
                if not quiet:
                    print "  * adding lookat movement to previous action - deleting old positions"
                addupStorage = [[0,0,0] for i in range(0,nsteps)]
            elif applymode == "addup":
                if not quiet:
                    print "  * adding lookat movement to previous action - adding up to old positions"
        else:
            addupStorage = [[0,0,0] for i in range(0,nsteps)]

        nBefore = len(self.lookat)

        if self.lastLookat() != start and applymode != "addup":
            self.lookat.append(start)
            nsteps = nsteps-1

        appendArray = []
        # ---------------------------------------------------------------------------------------------
        if stepmode=="linear":
            stepsizes = [(end[i]-start[i])*1./nsteps for i in range(0,3)]
            appendArray = [[start[i] + j*stepsizes[i] for i in range(0,3)] for j in range(1,nsteps+1)]
        # ---------------------------------------------------------------------------------------------
        elif stepmode=="log" or stepmode=="logrev":
            # some dirty tweaking here to have something well behaved (don't bother figuring out why and what :p)
            stepsizes = [math.pow(10,(math.log10(0.001)-math.log10(1.001))*1./nsteps)
                         if (start[i] != end[i]) else 1 for i in range(0,3)]
            if stepmode=="log":
                steplist = [[1-(stepsizes[i]**j - 0.001) if (start[i] != end[i]) else 1 for i in range(0,3)] 
                            for j in range(1,nsteps+1)]
            elif stepmode=="logrev":
                steplist = [[stepsizes[i]**j - 0.001 if (start[i] != end[i]) else 1 for i in range(0,3)] 
                            for j in range(0,nsteps)]
                steplist.reverse()
            steplist[-1] = [1.,1.,1.]
            appendArray = [[(end[i]-start[i])*steplist[j][i] for i in range(0,3)] 
                           for j in range(0,nsteps)]
        # ---------------------------------------------------------------------------------------------
        elif stepmode=="sin" or stepmode=="sinrev":
            amplitudes = [(end[i]-start[i])*1. for i in range(0,3)]
            stepsize = math.pi/(2*nsteps)
            xarray = [j*stepsize for j in range(1,nsteps+1)]
            if stepmode=="sin":
                yarray = [math.sin(x) for x in xarray]
            elif stepmode=="sinrev":
                yarray = [1-math.cos(x) for x in xarray]
            for j in range(0,len(yarray)):
                appendArray.append([start[i] + yarray[j]*amplitudes[i] for i in range(0,3)])          
        # ---------------------------------------------------------------------------------------------
        else:
            print "# ERROR : mode not recognized ["+mode+"]"
            exit(-1)
        # ---------------------------------------------------------------------------------------------

        if applymode=="addup":
            appendArray = [[addupStorage[j][i] + appendArray[j][i] for i in range(0,3)] for j in range(0,nsteps)]

        self.lookat.extend(appendArray)

        self.lastAction = dict()
        self.lastAction["action"] = "addLookatMovement"
        self.lastAction["nBefore"] = nBefore
        self.lastAction["nAfter"] = len(self.lookat)
        self.lastAction["expandedLists"] = []
        self.lastAction["modifiedLists"] = []
        if applymode == "expand":
            self.lastAction["expandedLists"] = ["lookat"]
        else:
            self.lastAction["modifiedLists"] = ["lookat"]

    def addSkyMovement(self, start, angle, nsteps, axis="dummy",applymode="expand"):
        """ 
        basic routine to add change in the sky vector
        - stepmode = 'linear' or 'log'
        - applymode = 'expand' or 'modify'
          'expand' just adds new sky vectors
          'modify' applies sky movement to previous action (if exists), possibly overwriting previously 
          generated sky vectors (so input variable nsteps is not used)
        """

        angle = angle*math.pi/180

        if applymode=="modify" and self.lastAction!=False:
            print "  * adding sky vector movement to previous action"
            index = self.lastAction["nBefore"]
            self.eraseList(index, self.sky)  
            nsteps = self.lastAction["nAfter"] - self.lastAction["nBefore"]

        if self.lastSky() != start:
            self.sky.append(start)
            nsteps = nsteps-1

        # testversion, just to quickly get some animation :p
        stepsize = angle*1./nsteps
        angle = 0
#        for j in range(0,nsteps):
#            angle += stepsize
#            newX = math.cos(angle)
#            newY = math.sin(angle)
#            self.sky.append([newX, newY,0])
        for j in range(0,nsteps):
            angle += stepsize
            newX = -math.sin(angle)
            newZ = math.cos(angle)
            self.sky.append([newX, 0, newZ])

    def addFidxProgress(self, start, end, nsteps, applymode="expand"):
        """ 
        basic routine to add progress in the file index (moving through time in sim)
        - applymode = 'expand' or 'modify'
          'expand' just adds new file indices
          'modify' applies file index progress to previous action (if exists), possibly overwriting previously 
          generated file indices (so input variable nsteps is not used)
        """

        if applymode=="modify" and self.lastAction!=False:
            print "  * adding fidx progress to previous action"
            index = self.lastAction["nBefore"]
            self.eraseList(index, self.fidx)  
            nsteps = self.lastAction["nAfter"] - self.lastAction["nBefore"]

        if self.lastFidx() != start:
            self.fidx.append(start)
            nsteps = nsteps-1

        stepsize = (end - start)*1./nsteps

        step = 0
        for j in range(0,nsteps):
            step += stepsize
            self.fidx.append(start+int(step))

    def addCamposOffset(self, offset, quiet=False):
        self.addCamposMovement(start=offset, end=offset, nsteps=0, stepmode="linear", 
                               applymode="addup", quiet=quiet)

    def padCampos(self, nsteps=None):
        if nsteps==None:
            nsteps = self.maxListLength() - len(self.campos)
        for j in range(0,nsteps):
            self.campos.append(self.lastCampos())

    def padLookat(self, nsteps=None):
        if nsteps==None:
            nsteps = self.maxListLength() - len(self.lookat)
        for j in range(0,nsteps):
            self.lookat.append(self.lastLookat())
    
    def padSky(self, nsteps=None):
        if nsteps==None:
            nsteps = self.maxListLength() - len(self.sky)
        for j in range(0,nsteps):
            self.sky.append(self.lastSky())

    def padFidx(self, nsteps=None):
        if nsteps==None:
            nsteps = self.maxListLength() - len(self.fidx)
        for j in range(0,nsteps):
            self.fidx.append(self.lastFidx())

    def maxListLength(self):
        return max(len(self.campos), len(self.lookat), len(self.sky), len(self.fidx))

    def maxValues(self, targetList):
        try:
            dim = len(targetList[0])
            return max([max(item) for item in targetList])
        except:
            try :
                targetList[0]
                return max(targetList)
            except:
                print "# ERROR : maxValues : targetList is empty"
                exit(-1)

    def eraseList(self, index, targetList):
        if index >= len(targetList):
            print "# WARNING : erasing : index",index,"larger than target list's length ("+str(len(targetList))+")"
            return
        deletedSlice = targetList[index:]
        del targetList[index:]
        return deletedSlice

    def padLists(self, nsteps=None):
        self.padCampos(nsteps=nsteps)
        self.padLookat(nsteps=nsteps)
        self.padSky(nsteps=nsteps)
        self.padFidx(nsteps=nsteps)

    def plotPath(self, maxCircle=5, plotLookat=False, plotLineOfSight=False, plotSky=False, 
                 plotSkyEvery=None, plotLosEvery=None):
        """ routine to plot the path in 3D """

        # setting up matplotlib stuff
        import matplotlib.pyplot as plt
        import mpl_toolkits.mplot3d.axes3d as plt3d
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax3d = plt3d.Axes3D(fig)

        # plotting campos data
        X = [campos[0] for campos in self.campos]
        Y = [campos[1] for campos in self.campos]
        Z = [campos[2] for campos in self.campos]
        ax3d.scatter(X,Y,Z, s=100, marker='d', color='red')
        ax3d.plot3D(X,Y,Z)

        # plotting lookat data
        if plotLosEvery == None:
            plotLosEvery = int((len(self.campos))*1./25)
        if plotLookat:
            lookatX = [lookat[0] for lookat in self.lookat]
            lookatY = [lookat[1] for lookat in self.lookat]
            lookatZ = [lookat[2] for lookat in self.lookat]
            ax3d.scatter(lookatX,lookatY,lookatZ, s=100, marker='d', color='red')
        if plotLineOfSight:
            losX = [[x,lookatx] for x,lookatx in zip(X,lookatX)]
            losY = [[y,lookaty] for y,lookaty in zip(Y,lookatY)]
            losZ = [[z,lookatz] for z,lookatz in zip(Z,lookatZ)]
            counter=0
            oldfidx=""
            for x,y,z, fidx in zip(losX,losY,losZ,self.fidx):
                if fidx!=oldfidx:
                    newFidx=True
                if (int(counter*1./plotLosEvery)-counter*1./plotLosEvery == 0.0 or newFidx==True):
                    if newFidx:
                        color="yellow"
                        linewidth="10"
                    else:
                        color="red"
                        linewidth=2
                    ax3d.plot3D(x,y,z,color=color,linewidth=linewidth)
                    newFidx=False
                oldfidx=fidx
                counter += 1

        # plotting sky vectors
        if plotSkyEvery == None:
            plotSkyEvery = int((len(self.campos))*1./25)
        if plotSky:
            vectorsize = self.maxValues(self.campos)/10.
            skyX = [sky[0]*vectorsize for sky in self.sky]
            skyY = [sky[1]*vectorsize for sky in self.sky]
            skyZ = [sky[2]*vectorsize for sky in self.sky]
            vectorx = [[x,x+skyx] for x,skyx in zip(X,skyX)]
            vectory = [[y,y+skyy] for y,skyy in zip(Y,skyY)]
            vectorz = [[z,z+skyz] for z,skyz in zip(Z,skyZ)]
            counter=0
            for x,y,z in zip(vectorx,vectory,vectorz):
                if (int(counter*1./plotSkyEvery)-counter*1./plotSkyEvery == 0.0):
                    ax3d.plot3D(x,y,z,color="green")
                counter += 1
            

        # setting plot properties/additions
        limit = int(2*max(maxCircle, self.maxValues(self.campos)))
        # axes
        Range = range(-limit,limit)
        zero = [0 for x in Range]
        ax3d.plot3D(Range,zero,zero,c='0.5')
        ax3d.plot3D(zero,Range,zero,c='0.5')
        ax3d.plot3D(zero,zero,Range,c='0.5')
        # circles
        for rad in range(1,maxCircle):
            angles = range(0,360)
            X = [rad*math.cos(angle*math.pi/180) for angle in angles]
            Y = [rad*math.sin(angle*math.pi/180) for angle in angles]
            Z = [0 for angle in angles]
            ax3d.plot3D(X,Y,Z,c='0.85')
        # limits
        ax3d.set_xlim3d(-limit,limit)
        ax3d.set_ylim3d(-limit,limit)
        ax3d.set_zlim3d(-limit,limit)
        
        plt.show()

    def writePath(self):

        #print len(self.campos), len(self.lookat), len(self.sky), len(self.fidx)
        self.padLists()
        sys.stdout.flush()

        output = open("geometry_"+self.filename,'w')
        print "* writing path to","geometry_"+self.filename

        for j in range(0,len(self.campos)):
            #fidx = str(self.fidx[j])
            #if self.fidx[j] < 100:
            #    fidx = "0" + fidx
            #if self.fidx[j] < 10:
            #    fidx = "0" + fidx
            output.write(str(self.campos[j][0])+" "+ str(self.campos[j][1])+" "+ str(self.campos[j][2])+" "+
                         str(self.lookat[j][0])+" "+ str(self.lookat[j][1])+" "+ str(self.lookat[j][2])+" "+
                         str(self.sky[j][0])+" "+ str(self.sky[j][1])+" "+ str(self.sky[j][2])+" "+
                         "%04d"%self.fidx[j] + "\n")
    
        output.close()






if __name__ == "__main__":

	filename = "sim62002" # fly around gas 3

	#path = pathGenerator(filename = filename, campos = [0,7.5,0], lookat = [0,0,0], sky = [0,0,1], fidx = 35)
	path = pathGenerator(filename = filename, campos = [0,5,0], lookat = [0,0,0], sky = [0,0,1], fidx = 0)
	#path.addSimpleLineSegment(startCampos="continue", endCampos=[0,10,0], nsteps=100,mode="linear")
   	path.addFidxProgress(start=0, end=130, nsteps=129, applymode="expand")   
        #path.addHelix(radius=5, axis="z", rotation=2, direction="righthanded", nsteps=200, displacement=0)
        #path.addHelix(radius=15, axis="z", rotation=2, direction="righthanded", nsteps=200, displacement=0, center=[-0.025,-0.7,0.6])
   	#path.addFidxProgress(start=9, end=35, nsteps=24, applymode="expand")   
#    path.addFidxProgress(start=0, end=16, nsteps=2100, applymode="modify")
        #path.addSimpleLineSegment(startCampos="continue", endCampos=[0,5,0], nsteps=100,mode="linear")
	#path = pathGenerator(filename = filename, campos = [0,5,0], lookat = [0,0,0], sky = [0,0,1], fidx = 478)
        #path.addHelix(radius=5, axis="z", rotation=1, direction="righthanded", nsteps=100, displacement=0, center=[-0.025,-0.7,0.6])
   	#path.addFidxProgress(start=478, end=559, nsteps=86, applymode="expand")   
	#path.addSimpleLineSegment(startCampos="continue", endCampos=[1,0,-12], nsteps=250,mode="sin")
#    path.addSimpleLineSegment(startCampos="continue", endCampos=[-1,1,0], nsteps=250,mode="sinrev")
#    path.addSimpleLineSegment(startCampos="continue", endCampos=[6,0,0], nsteps=250,mode="sin")
#    path.addHelix(radius=6, axis="z", rotation=5, direction="righthanded", nsteps=2500, displacement=0, center=[0,0,0])
##    path.addFidxProgress(start=16, end=36, nsteps=2500, applymode="modify")
#    path.addSimpleLineSegment(startCampos="continue", endCampos=[4,0,0], nsteps=100,mode="sin")
#    path.addHelix(radius=4, axis="z", rotation=5, direction="righthanded", nsteps=2500, displacement=0, center=[0,0,0])
#    path.addFidxProgress(start=36, end=56, nsteps=2500, applymode="modify")

     	path.writePath()
	#path.plotPath(plotLookat=True,plotLineOfSight=False,plotSky=False)



