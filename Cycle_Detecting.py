import sys
import numpy as np 
import os
from  Detect_C import Simple_INFO



class Ring_Detecting(Simple_INFO.Single_Bond):  
    def __init__(self,NAtom,BO,IAn,ICoor):
        Simple_INFO.Single_Bond.__init__(self,NAtom,BO,IAn,ICoor)
#        self.NAtom        = NAtom
#        self.BO           = BO
#        self.IAn          = IAn
#        self.Terminal_Atom= [1,9,17,35,53]
#        self.Halogen_Atom = [9,17,35,53]


    def detect_cycle(self):
        layer_MAX = 15
        layer = 1
        I = 0
        J = 0 
        First_Atom = 0  # The first atom of a cycle 
        I_Cycle = []  # the inital cycle is zero
        Initial = []
        Cycle_ONE = []
        Final_Cycle = []
        Ring_TMP = []
        while layer >= 1 :
            while I < self.NAtom:
                if ( First_Atom != I ):
#                    Initial = []   # if the initial atom has changed, the initial should be empty
                    First_Atom = I      # the initial initial atom changes 
                if (layer == 1) :
                    Cycle_ONE.append(I)
                    if (layer != len(Cycle_ONE)):
                        Cycle_ONE = [I]
                if (self.IAn[Cycle_ONE[-1]] in self.Terminal_Atom):  #H Atom and halogen Atoms can't be the skeleton of ring
                    I = I + 1
                    continue 
                else: 
                    while J < self.NAtom: 
                        if (layer > layer_MAX):
                            layer = layer - 1
                            J=Cycle_ONE[-1] + 1
                            Cycle_ONE.pop(-1)
                            continue
                        if (layer == 1):
                            replace_1 = 2
                        else:
                            replace_1 = layer
                        if (self.IAn[J] in self.Terminal_Atom ): #H Atom and halogen Atoms can't be the skeleton of ring
                            J = J + 1
                            continue
                        if ((self.BO[Cycle_ONE[-1]][J] != 0) and (J != Cycle_ONE[replace_1-2]) ):
                            if (layer > 3) :
                                if (J in Cycle_ONE[1:]) :
                                    J =J + 1
                                    continue
                            if (J == Cycle_ONE[0]):
                                TMP = Cycle_ONE[:]
                                I_Cycle.append(TMP)
    #                            MAT = I_Cycle
                                J = J + 1
                                continue
                            else:
                                layer =layer + 1
                                Cycle_ONE.append(J)
                                J = 0
                                continue
                        else:
                            J = J + 1
                            continue
                    layer = layer - 1 
                    if (layer == 0):
                        I = I + 1
                        J = 0
#                        I_Cycle.append(Initial)   #
                        layer = layer + 1
#                        Initial = []
                    else: 
                        J = Cycle_ONE[-1]+1
                    Cycle_ONE.pop(-1)   # the connected atom shouldn't be H Atom or halogen Atoms 
            layer = layer -1   
#*******************The following code is aim to sift the repeat rings **************
        for Ring in I_Cycle:
            Ring_sort = sorted(Ring)
            if (Ring_sort not in Ring_TMP):
                Final_Cycle.append(Ring)
                Ring_TMP.append(Ring_sort)    
        
#***************************************************************************************
        M = 0
        N = 0
        J = 0
        List_Atom_in_Circle_repeat_Value = [] 
        List_Circle_Atom = []
        while M < len(Final_Cycle):
            while N < len(Final_Cycle[M]):
                List_Atom_in_Circle_repeat_Value.append(Final_Cycle[M][N])
                N = N + 1
            N = 0
            M = M + 1

        List_Atom_in_Circle = list(set(List_Atom_in_Circle_repeat_Value))
        while  J < self.NAtom:
            if J in List_Atom_in_Circle:
                List_Circle_Atom.append(1)
            else:
                List_Circle_Atom.append(0)
            J = J + 1


  #      return I_Cycle,Final_Cycle,List_Atom_in_Circle,List_Circle_Atom
  #      return I_Cycle,Final_Cycle
        return Final_Cycle, List_Circle_Atom

if __name__ == "__main__" :
    if (len(sys.argv) < 2) :
        print ("Usage: python cycle.py BO IAn")
        os._exit()
    
    Path1 = sys.argv[1]
    Path2 = sys.argv[2]
    BO = np.loadtxt(Path1)
#    IAn = np.loadtxt(Path2)
    IAn = (np.loadtxt(Path2))[:,0]
    ICoor = (np.loadtxt(Path2))[:,1:]
    NAtom = len(IAn)
    X = Ring_Detecting(NAtom,BO,IAn)
    _I_Cycle, _Final_Cycle,_List_Atom_in_Circle,_List_Circle_Atom = X.detect_cycle()
    print (_I_Cycle)
    print(_Final_Cycle)
    print(_List_Atom_in_Circle)
    print(_List_Circle_Atom)


    







