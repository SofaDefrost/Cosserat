import Sofa
from splib3.numerics import Quat


class FingerController(Sofa.Core.Controller):
    """
        Implements the AnimationManager as a PythonScriptController
    """

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.rigidBaseMO = args[0]
        self.rateAngularDeformMO = args[1]

        self.rate = 0.2
        self.angularRate = 0.02
        return

    def _extracted_from_onKeypressedEvent_10(self, qOld, posA, angularRate):
        qNew = Quat.createFromEuler([0., angularRate, 0.], 'ryxz')
        qNew.normalize()
        qNew.rotateFromQuat(qOld)
        for i in range(4):
            posA[0][i+3] = qNew[i]

    def onKeypressedEvent(self, event):
        key = event['key']
        if ord(key) == 19:  # up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(4):
                    qOld[i] = posA[0][i+3]

                self._extracted_from_onKeypressedEvent_10(
                    self, qOld, posA, self.angularRate)

        if ord(key) == 21:  # down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(4):
                    qOld[i] = posA[0][i+3]

                self._extracted_from_onKeypressedEvent_10(
                    qOld, posA, -self.angularRate)

        # Pull the cable base
        if ord(key) == 18:  # left
            print("left: the cable is pulled")
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] -= self.rate
                print(f' The new position is : {posA[0]}')
        # release the cable base
        if ord(key) == 20:  # right
            print("left: the cable is released")
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] += self.rate
                print(f' The new position is : {posA[0]}')