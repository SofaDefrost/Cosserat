#!/usr/bin/env python
# -*- coding: utf-8 -*-

import Sofa

class controller(Sofa.PythonScriptController):

    def initGraph(self, node):
            self.node = node

    def onKeyPressed(self,c):
            inputvalue = self.node.getObject('RigidBaseMO').findData('resPosition')

            if (c == "+"):
                position = inputvalue.resPosition[0]
                displacement = inputvalue.resPosition[0][0] + 0.1
                print("++++++++++++++++++++++++++++++++++++++++++++++++++>"+str(displacement) +" "+str(position))
                #inputvalue.resPosition = str(displacement)

            elif (c == "-"):
                position = inputvalue.resPosition[0]
                displacement = inputvalue.resPosition[0][0] - 0.1
                print("-------------------------------------------------------> "+str(displacement) +" position :"+str(position))
                if(displacement < 0):
		  displacement = 0
               #inputvalue.value = str(displacement)
