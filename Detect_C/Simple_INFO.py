import numpy as np 
import sys
import os

class Single_Bond:
    def __init__(self,NAtom,BO,IAn,ICoor): # Nomal_Bond is the number of bond that the initial atom aonnected with other, for instance, C is 4,O is 2
        self.NAtom        = NAtom
        self.BO           = BO
        self.IAn          = IAn
        self.ICoor        = ICoor
        self.Terminal_Atom= [1,9,17,35,53]  #Terminal atoms are the atom which attatch with one atom with only one single bond, like H, F, Cl, Br, I 
        self.Halogen_Atom = [9,17,35,53]    
#       **********************
        self.C_Nuclei_Num = 6
        self.H_Nuclei_Num = 1
        self.O_Nuclei_Num = 8
        self.N_Nuclei_Num = 7
        self.B_Nuclei_Num = 5 
        self.S_Nuclei_Num = 16
        self.P_Nuclei_Num = 15
        self.F_Nuclei_Num = 9
        self.Cl_Nuclei_Num = 17
        self.Br_Nuclei_Num = 35
        self.I_Nuclei_Num = 53
        self.As_Nuclei_Num = 33
        self.Sb_Nuclei_Num = 51
#       **********************
        self.C_Normal_Bond = 4
        self.Atom_Has_Been_Count = []
        self.Dis = self.Distance()
        self.BO  = self.BO_Self_Correction()
        self.List_Heavy_Atom = self.Heavy_Atom()


    def Attach(self):   # this function is aim to detect the connecting situation of each atom 
        I = 0
        J = 0
        Num_Attached_Atom = 0
        Num_Attached_Bond = 0 
        List_Num_Attached_Atom = []
        List_Num_Attached_Bond = []
        while I < self.NAtom :
            while J < self.NAtom:
                Num_Attached_Bond = Num_Attached_Bond + self.BO[I][J]  # add all the Bond order together 
                if (self.BO[I][J] !=0):
                    Num_Attached_Atom = Num_Attached_Atom + 1   # detect how many atom the Initial atom attachs with
                J = J + 1
            List_Num_Attached_Atom.append(int(Num_Attached_Atom))
            List_Num_Attached_Bond.append(int(Num_Attached_Bond))
            Num_Attached_Atom = 0
            Num_Attached_Bond = 0
            I = I + 1
            J = 0 
        return List_Num_Attached_Atom,List_Num_Attached_Bond
           

    def Hybrid(self,Attached_Bond,Attached_Atom,List_Surrounding_Atom):
        I = 0 
        J = 0
        List_Hybrid = []
        while I < self.NAtom:
#********************** This part code is for detecting the hybird method of C
            if (self.IAn[I] == self.C_Nuclei_Num):
                if (Attached_Bond[I] == 4 and Attached_Atom[I] == 4):
                    List_Hybrid.append(3)
                    I = I + 1
                    continue
                if (Attached_Bond[I] == 4 and Attached_Atom[I] == 3):
                    List_Hybrid.append(2)
                    I = I + 1
                    continue
                if (Attached_Bond[I] == 4 and Attached_Atom[I] == 2):
                    List_Hybrid.append(1)
                    I = I + 1
                    continue
                if (Attached_Bond[I] >= 5 and Attached_Atom[I] == 3): # Attached Bond is 4 under the normal situation, however, due to the special structure of benzene whose C atom is seem as connecting with another C atom with Double bond , which means that the Attached Bond exceeds 4, equals 5
                    List_Hybrid.append(2)
                    I = I + 1
                    continue
                if (Attached_Bond[I] == 3 and Attached_Atom[I] == 3):
                    _Dihedral_Angle = self.Dihedral_Angle(List_Surrounding_Atom[I])
#                    print(I,_Dihedral_Angle,"(((((((((((((((((((")
                    if (_Dihedral_Angle >= 165 and _Dihedral_Angle <= 180):
                        List_Hybrid.append(2)
                        I = I + 1
                        continue
                    
#*********************CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC**************************
#################This part code is aim to anaysis the hybird method of N and Atom
            if (self.IAn[I] == self.N_Nuclei_Num or self.IAn[I] == self.P_Nuclei_Num):
                if (Attached_Bond[I] == 3 and Attached_Atom[I] == 3):
                    _Dihedral_Angle = self.Dihedral_Angle(List_Surrounding_Atom[I])
#                    print(_Dihedral_Angle,I,"NNNNNNNNNNNNNNNNNNNNNN")
                    if (_Dihedral_Angle >= 165 and _Dihedral_Angle <= 180):
                        List_Hybrid.append(2)
                    else :
                        List_Hybrid.append(3)
                    I = I + 1
                    continue
                if (Attached_Bond[I] == 3 and Attached_Atom[I] == 2):
                    List_Hybrid.append(2)
                    I = I + 1
                    continue
                if (Attached_Bond[I] == 4 and Attached_Atom[I] == 4):
                    List_Hybrid.append(3)
                    I = I + 1
                    continue
                if (Attached_Bond[I] > 3 and Attached_Atom[I] == 3):
                    List_Hybrid.append(2)
                    I = I + 1
                    continue
                if (Attached_Bond[I] == 3 and Attached_Atom[I] == 1):
                    List_Hybrid.append(1)
                    I = I + 1
                    continue
#***********************NNNNNNNNNNNNNNNNNNNNNNNNNNN*****************************
#********************* This part code is in order to analysis the hybird method of O ******

            if (self.IAn[I] == self.O_Nuclei_Num or self.IAn[I] == self.S_Nuclei_Num):
                if (Attached_Bond[I] == 2 and Attached_Atom[I] == 2):
                    List_Hybrid.append(3)
                    I = I + 1
                    continue
                if (Attached_Bond[I] == 2 and Attached_Atom[I] == 1):
                    List_Hybrid.append(2)
                    I = I + 1
                    continue
                if (Attached_Bond[I] == 1 and Attached_Atom[I] == 1):
                    List_Hybrid.append(3)
                    I = I + 1
                    continue
                if (Attached_Bond[I] == 3 and Attached_Atom[I] == 3):   
                    List_Hybrid.append(3)
                    I = I + 1
                    continue
#################OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO***********************************
            I = I + 1
            List_Hybrid.append(0)
        return List_Hybrid


    def Surrounding_Atom(self):
        List_Surrounding_Atom = []
        I = 0
        J = 0
        TMP = []
        while I < self.NAtom:
            TMP.append(I)
            while J < self.NAtom:
                if (self.BO[I,J] != 0 ):
                    TMP.append(J)
                J = J + 1
            List_Surrounding_Atom.append(TMP[:])
            TMP = []
            J = 0
            I = I + 1
        return  List_Surrounding_Atom



    def Num_of_terminal_atom_Surrounding(self):   # this method is designed for analysis the number of terminal_atom (H, F, Cl, Br) surrounding the central atom like -CF3, so the number is 3 
        List_Num_of_terminal_atom_Surrounding = []
        I = 0
        J = 0
        Icount = 0 
        while I < self.NAtom:
            while J < self.NAtom:
                if (self.BO[I][J] == 1 and self.IAn[J] in self.Terminal_Atom ):
                    Icount = Icount + 1
                J = J + 1
            List_Num_of_terminal_atom_Surrounding.append(Icount)
            I = I + 1
            Icount = 0 
            J = 0
        return List_Num_of_terminal_atom_Surrounding

    


    def Dihedral_Angle(self,List_Surrounding_Atom): 
        # this method is based on the formula in the wiki https://en.wikipedia.org/wiki/Dihedral_angle
        p2 = self.ICoor[List_Surrounding_Atom[0]]
        p1 = self.ICoor[List_Surrounding_Atom[1]]
        p3 = self.ICoor[List_Surrounding_Atom[2]]
        p4 = self.ICoor[List_Surrounding_Atom[3]]

        b1 = -1.0*(p1 - p2)
        b2 = p3 - p2
        b3 = p4 - p3

        b1xb2 = np.cross(b1,b2)
        b2xb3 = np.cross(b2,b3)

        b1xb2_x_b2xb3 = np.cross(b1xb2, b2xb3)

        y = np.dot(b1xb2_x_b2xb3, b2)*(1.0/np.linalg.norm(b2))
        x = np.dot(b1xb2, b2xb3)
        return abs(np.degrees(np.arctan2(y,x)))




    def BO_Self_Correction(self): # this method is desiged for correcting the BO matrix, since there are some mistake in BO matrix, for instance, when Br Atom connects with C atom, BO matrix demonstrates that C and Br Atom dont connect with each other. so we should correct the BO matrix and return a realtively precise BO matrix.
        I = 0
        J = 0
        K = 0
        BO = self.BO
        if (self.Br_Nuclei_Num in self.IAn or self.I_Nuclei_Num in self.IAn):
            while I < self.NAtom:
                if (self.IAn[I] == self.Br_Nuclei_Num ):
#                    print(I,"XXXXXXXXXXXXXXXXXXXXXX")
                    _The_Closest_Atom = self.Finding_The_Closest_Atom(I,self.Dis[I])  # find the atom which is closest to the Br atom
                    BO[I][_The_Closest_Atom] = 1
                    BO[_The_Closest_Atom][I] = 1
                if (self.IAn[I] == self.I_Nuclei_Num ):
#                    print(I,"XXXXXXXXXXXXXXXXXXXXXX")
                    _The_Closest_Atom = self.Finding_The_Closest_Atom(I,self.Dis[I])  # find the atom which is closest to the Br atom
                    BO[I][_The_Closest_Atom] = 1
                    BO[_The_Closest_Atom][I] = 1
                I = I + 1
        while J < self.NAtom:
            while K < self.NAtom:
                if (BO[J][K] == 5 ): # this is Hydrogen bond
                    BO[J][K] = 0
                K = K + 1
            K = 0
            J = J + 1
        return BO    
           



    def Distance(self): # this function is aimed at calculating the distance between two atom 
        Dis = []
        I  = 0
        J  = 0 
        while I < self.NAtom:
            Dis_Temp = []
            while J < self.NAtom:
                Dis_Atom = np.sqrt(np.sum(np.square(self.ICoor[I] - self.ICoor[J])))
                Dis_Temp.append(Dis_Atom)
                J = J + 1
            J = 0
            I = I + 1
            Dis.append(Dis_Temp)
        return Dis



    def Finding_The_Closest_Atom(self,Central_Atom_Num,List_Distance): # this function is aimed at finding the atom which is closest to the Br atom(except H Atom)
        I = 0 
        List_Candidate_Atom = []
        List_Distance_Between_Candidate_Atom_And_Central_Atom = []
        while I < self.NAtom:
            if ( I != Central_Atom_Num and self.IAn[I] != self.H_Nuclei_Num):
                List_Candidate_Atom.append(I)
                List_Distance_Between_Candidate_Atom_And_Central_Atom.append(self.Dis[Central_Atom_Num][I])
            I = I + 1
        The_Closet_Atom_Num = List_Candidate_Atom[List_Distance_Between_Candidate_Atom_And_Central_Atom.index(min(List_Distance_Between_Candidate_Atom_And_Central_Atom))]
        return The_Closet_Atom_Num
                
        
    

    def Heavy_Atom(self):  # this function is trying to find all heavy atom (Non-hydrogen Atom) 
        List_Non_Hydrogen_Atom = []
        I = 0
        while I < self.NAtom:
            if (self.IAn[I] != 1):
                List_Non_Hydrogen_Atom.append(I)
            I = I + 1
        return List_Non_Hydrogen_Atom
        



if __name__ == "__main__":
    if (len(sys.argv) < 2) :
        print ("Usage: python cycle.py BO IAn")
        os._exit()

    Path1 = sys.argv[1]
    Path2 = sys.argv[2]
    BO = np.loadtxt(Path1)
    print(BO[18])
    IAn = (np.loadtxt(Path2))[:,0]
    ICoor = (np.loadtxt(Path2))[:,1:]
    NAtom = len(IAn)
    X =  Single_Bond(NAtom,BO,IAn,ICoor)
    List_Num_Attached_Atom, List_Num_Attached_Bond = X.Attach()
    List_Surrounding_Atom                              = X.Surrounding_Atom()
    List_Hybrid                                        = X.Hybrid(List_Num_Attached_Bond,List_Num_Attached_Atom,List_Surrounding_Atom)
    List_Num_of_terminal_atom_Surrounding              = X.Num_of_terminal_atom_Surrounding()
#    print(List_Num_Attached_Atom, List_Num_Attached_Bond)
    print (List_Hybrid[12])
    print (List_Hybrid[10])
    print (List_Hybrid[18])
    print (List_Hybrid[27])
    print(X.BO[12])
