import Sofa
import numpy as np
from splib3.constants import Key
import time
import csv


class DirectController(Sofa.Core.Controller):

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.name = "DirectController"
        self.cableConstraint =  kwargs['cableConstraint']
        self.constantForceField = kwargs['constantForceField']
        self.rootNode =  kwargs['rootNode']
        self.cfmax = kwargs['cfMax']
        self.pfmax = kwargs['pfMax']
        self.risingTime = kwargs['risingTime']
        self.time = 0.0
        self.dt = self.rootNode.dt.value
        self.cableForce = 0.0
        self.ponctualForce = 0.0

        self.saveFlag = False
        self.saveTime = 30.0

        self.mechanicalObject = kwargs['mechanicalObject']
        self.prevMesTime = time.time()
        self.listData = []

    def onAnimateEndEvent(self,event):

        self.time += self.dt

        if self.time < self.saveTime:

            if self.time < self.risingTime:
                self.cableForce += self.cfmax*self.dt/self.risingTime
                self.ponctualForce += self.pfmax*self.dt/self.risingTime
            else:
                self.cableForce = self.cfmax
                self.ponctualForce = self.pfmax

            print(f'The controller is called with force : {self.constantForceField.totalForce}')
            self.cableConstraint.value = [self.cableForce]
            # self.constantForceField.totalForce.value = [0.0,0.0,0.0,0.0,-self.ponctualForce,0.0]
            # self.constantForceField.totalForce.value = [0.0,-self.ponctualForce,0.0,0.0,0.0,0.0]



            posEff = self.mechanicalObject.position.value[-1]
            elapsedTime = time.time()-self.prevMesTime
            self.prevMesTime = time.time()
            #print("ElapsedTime: "+str(elapsedTime))
            #print("Rate: "+str(1.0/elapsedTime))

            self.listData.append([self.time,elapsedTime] + [posEff[0],posEff[1],posEff[2]])

        else:

            if not self.saveFlag:
                print(self.listData)
                # with open('fullCosserat_Nfr2_Ns_17.csv', 'w', newline='') as f:
                # with open('fullFEM_large_rigidDisks_hexa_Nx_801.csv', 'w', newline='') as f:
                with open('FEMProj_large_rigidDisks_lc_2_Nfr_2_Ns_2_noVisu.csv', 'w', newline='') as f:
                    # using csv.writer method from CSV package
                    write = csv.writer(f, delimiter=',',
                                        quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    for k in range(0, len(self.listData)):
                        write.writerow(self.listData[k])
                self.saveFlag = True
