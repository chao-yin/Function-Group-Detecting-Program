#from Cycle_Detecting import Ring_Detecting
from Detect_C import Group_Detect
import sys
import numpy as np
import os



#class Mul_Ring(Ring_Detecting,Group_Detect.Mutiple_Bond_Functional_Group):
class Mul_Ring(Group_Detect.Mutiple_Bond_Functional_Group):
    def __init__(self,NAtom,BO,IAn,ICoor):
#        Ring_Detecting.__init__(self,NAtom,BO,IAn)
        Group_Detect.Mutiple_Bond_Functional_Group.__init__(self,NAtom,BO,IAn,ICoor)
#        _,self.Ring_Set = Ring_Detecting.detect_cycle(self)

 
    def Resemble_No_Hetero_Ring_Group(self, Total_Num_Atom, Hybrid_Method, Majority_Atom_Type = 6): # this code can only detect Cyclopentyl(5,3),Cyclohexyl(6,3) and Phenyl(6,2)
        Icount = 0
        Num_Group = 0
        Set_Aimed_Group = []
        for Ring in self.Ring_Set:
            Num_Atom_of_Ring = len(Ring)
            if (Num_Atom_of_Ring == Total_Num_Atom):
#                print(Ring,"<<<<<<<<<<<<<<<<<<<<<<<<")
                for I in Ring:
#                    print(self.List_Hybrid[I],self.IAn[I])
                    if (self.List_Hybrid[I] == Hybrid_Method and self.IAn[I] == Majority_Atom_Type):
                        Icount = Icount + 1 
#                        print (I,"###################")
                     #   print(I,self.IAn[I])
#                print(Icount)
                if (Icount == Num_Atom_of_Ring):
                    
                    Num_Group = Num_Group + 1
                    for M in Ring:
                        self.Atom_Has_Been_Count.append(M)           #  adding the atom to the list of atom having been counted.
                    Set_Aimed_Group.append(Ring)
                Icount = 0
        return Set_Aimed_Group, Num_Group 
                    
    
    def Resemble_Single_Hetero_Ring_Group(self, Total_Num_Atom, Hybrid_Method, Hetero_Atom_Type, Num_Bond_Hetero_Atom_Attaching, Majority_Atom_Type = 6):
        Icount = 0 
        Num_Group = 0
        Num_Hetero_Atom = 0
        for Ring in self.Ring_Set:
            Num_Atom_of_Ring = len(Ring)
#            print("_____________________")
            if (Num_Atom_of_Ring == Total_Num_Atom):
#                print(Ring,"OOOOOOOOOOOOOOO")
                for I in Ring:
                    if (self.IAn[I] == Hetero_Atom_Type):
#                        print(I,self.IAn[I],"CCCCCCCCCCCC" ) 
                        Hetero_Atom = int(I)
                        Num_Hetero_Atom = Num_Hetero_Atom + 1
#                    print (self.List_Hybrid[I])
                    if (self.List_Hybrid[I] == Hybrid_Method and self.IAn[I] == self.C_Nuclei_Num):
#                        print(I,"????????????")
                        Icount = Icount + 1
#                print (Icount,"IIIIIIIIIIIIIII")
                if (Num_Hetero_Atom == 1):
#                    print (Hetero_Atom,int(self.IAn[Hetero_Atom]),Ring,"************")
#                    print(Icount)
                    if (Icount == Total_Num_Atom - 1):
#                        print(Ring,"!!!!!!!!!!")
                        if (self.List_Num_Attached_Atom[Hetero_Atom] == Num_Bond_Hetero_Atom_Attaching):
#                            print(Ring,"@@@@@@@@@@@@@@@@")
                            Num_Group = Num_Group + 1
                            for M in Ring:
                                self.Atom_Has_Been_Count.append(M)           #  adding the atom to the list of atom having been counted.
                Icount = 0
                Num_Hetero_Atom = 0
        return Num_Group

    
    


    def Resemble_Double_Hetero_Ring_Group(self, Membered_Rings, Central_Hetero_Atom_Type, Num_of_Atom_Central_Atom_Connected, Num_of_Bond_Central_Atom_Connected, Subsidiary_Hetero_Atom_Type, Num_of_Atom_Subsidiary_Atom_Connected, Num_of_Bond_Subsidiary_Atom_Connected, Hybrid_Method_Of_All_C, Distance_between_Central_And_Subsidiary_Atom, Central_And_Subsidiary_method = 0, Num_of_C_Abnormal_Hybridzation_Method = 0):  # Membered_Rings means which membered this ring is, Central_Hetero__Atom_Type means the type of central hetero atom, like the central hetero atom of po3 is 16(S), Subsidiary_Hetero_Atom_Type is another atom like po3 is 7(N), Hybrid_Method of_All_C is the hybrid method of all C, in po3 this value is 2 standing for sp2 hybridzation, Distance_between_Central_And_Subsidiary_Connecting means whether central hetero atom are connected with subsidary atom, 1 stands for conneting, 2 means that central and subsidiary atom are isolated by 1 atom.Central_And_Subsidiary_method stands for the conneted bond type, 0 stand for disconnected, 1 stands for single bond, default is 0. so the expression of po3 is (5,16,3,3,7,2,3,2,2), the number of arguments of this function is 11 containing 2 defaults. 
        Central_Hetero_Atom = []
        Subsidiary_Hetero_Atom = []
        C_Atom = []
        Num_of_targeted_Group = 0
        Set_target_Group = []
        for Ring in self.Ring_Set:
            Num_Atom_of_Ring = len(Ring)
            if ( Num_Atom_of_Ring == Membered_Rings ):
#                print(Ring,"AAAAAAAAAAAAAAAAA")
                for I in Ring:
                    if (self.IAn[I] == Central_Hetero_Atom_Type and self.List_Num_Attached_Atom[I] == Num_of_Atom_Central_Atom_Connected and self.List_Num_Attached_Bond[I] == Num_of_Bond_Central_Atom_Connected):   # search for the central atom which is S and connecting with 3 atoms and the connecting bond sum is 3 
                        Central_Hetero_Atom.append(I)
#                        print (I,"PPPPPPPPPPPPPPP", self.List_Num_Attached_Atom[I], self.List_Num_Attached_Bond[I])
                    if (int(self.IAn[I]) == Subsidiary_Hetero_Atom_Type and self.List_Num_Attached_Atom[I] == Num_of_Atom_Subsidiary_Atom_Connected and self.List_Num_Attached_Bond[I] == Num_of_Bond_Subsidiary_Atom_Connected ):  # search for the Subsidiary atom which is S and connecting with 3 atoms and the connecting bond sum is 3
                        Subsidiary_Hetero_Atom.append(I)
#                        print (I,"SSSSSSSSSSSSSSSSSS", self.List_Num_Attached_Atom[I], self.List_Num_Attached_Bond[I])
                    if (self.IAn[I] == self.C_Nuclei_Num and self.List_Hybrid[I] == Hybrid_Method_Of_All_C):
                        C_Atom.append(I) 
#                        print (I,"CCCCCCCCCCCCC")
                if (Central_Hetero_Atom_Type == Subsidiary_Hetero_Atom_Type and Num_of_Atom_Central_Atom_Connected == Num_of_Atom_Subsidiary_Atom_Connected and  Num_of_Bond_Central_Atom_Connected == Num_of_Bond_Subsidiary_Atom_Connected ):
                    Num_Central_Hetero_Atom = len(Central_Hetero_Atom)/2
                    Num_Subsidiary_Hetero_Atom = len(Subsidiary_Hetero_Atom)/2
                else:
                    Num_Central_Hetero_Atom = len(Central_Hetero_Atom)
                    Num_Subsidiary_Hetero_Atom = len(Subsidiary_Hetero_Atom)
                Num_C_Atom = len(C_Atom)

#                print(Num_Central_Hetero_Atom, Num_Subsidiary_Hetero_Atom,  Num_C_Atom,"*************", Membered_Rings-2-Num_of_C_Abnormal_Hybridzation_Method )
                if ( Num_Central_Hetero_Atom == 1 and Num_Subsidiary_Hetero_Atom == 1 and Num_C_Atom == Membered_Rings-2-Num_of_C_Abnormal_Hybridzation_Method ):
                    if (Central_Hetero_Atom_Type == Subsidiary_Hetero_Atom_Type and Num_of_Atom_Central_Atom_Connected == Num_of_Atom_Subsidiary_Atom_Connected and  Num_of_Bond_Central_Atom_Connected == Num_of_Bond_Subsidiary_Atom_Connected ):
                        Position_of_Central_Atom_in_rings = Ring.index(Central_Hetero_Atom[0])
                        Position_of_Subsidiary_Atom_in_rings = Ring.index(Subsidiary_Hetero_Atom[1])
                    else:
                        Position_of_Central_Atom_in_rings = Ring.index(Central_Hetero_Atom[0])
#                        print(Position_of_Central_Atom_in_rings,Central_Hetero_Atom[0])
                        Position_of_Subsidiary_Atom_in_rings = Ring.index(Subsidiary_Hetero_Atom[0])
#                        print(Position_of_Subsidiary_Atom_in_rings,Subsidiary_Hetero_Atom[0])
                    Distance_Central_and_Subsidiary = abs(Position_of_Central_Atom_in_rings-Position_of_Subsidiary_Atom_in_rings) 
#                    print("DDDDDDDDDD", Distance_Central_and_Subsidiary, Distance_between_Central_And_Subsidiary_Atom)
                    if (Distance_Central_and_Subsidiary > int(Membered_Rings/2)):
                        Distance_Central_and_Subsidiary = Membered_Rings - Distance_Central_and_Subsidiary    # since the ring is a circle, so the distance must be littler than half of the membered of ring. 
                    if ( Distance_Central_and_Subsidiary == Distance_between_Central_And_Subsidiary_Atom ):  # matching the distance between the central and Subsidiary atom
                        if ( Distance_between_Central_And_Subsidiary_Atom == 1 ):
                            if (self.BO[Central_Hetero_Atom[0]][Subsidiary_Hetero_Atom[0]] == Central_And_Subsidiary_method):  # if the central and Subsidiary are connecting with each other, then check the bond type of the bond
                                 Num_of_targeted_Group = Num_of_targeted_Group + 1
                                 Set_target_Group.append(Ring)
                                 for M in Ring:
                                     self.Atom_Has_Been_Count.append(M)           #  adding the atom to the list of atom having been counted.
                        else: 
                            Num_of_targeted_Group = Num_of_targeted_Group + 1 
                            Set_target_Group.append(Ring)
                            for M in Ring:
                                self.Atom_Has_Been_Count.append(M)           #  adding the atom to the list of atom having been counted.
                Central_Hetero_Atom = []
                Subsidiary_Hetero_Atom = []
                C_Atom = [] 		
        return Num_of_targeted_Group, Set_target_Group					





    def Five_Membered_Ring_With_Three_N_Atom(self):
        Central_Atom = []
        One_Subsidiary_Atom = []
        Another_Subsidiary_Atom = []
        Num_Target_Group = 0
        C_Atom = []
        for Ring in self.Ring_Set:
            Num_Atom_of_Ring = len(Ring)
            if (Num_Atom_of_Ring == 5 ):
       #         print(Ring)
                for I in Ring :
                    if (self.IAn[I] == self.N_Nuclei_Num and self.List_Num_Attached_Atom[I] == 3 and self.List_Num_Attached_Bond[I] == 4 ) :
                        Central_Atom.append(I)
#                        print (I,"Central")
                    if (self.IAn[I] == self.N_Nuclei_Num and self.List_Num_Attached_Atom[I] == 2 and self.List_Num_Attached_Bond[I] == 3 ):
                        One_Subsidiary_Atom.append(I)
#                        print (I,"One_Subsidiary")
                    if (self.IAn[I] == self.N_Nuclei_Num and self.List_Num_Attached_Atom[I] == 3 and self.List_Num_Attached_Bond[I] == 3 ):
#                    if (self.IAn[I] == self.N_Nuclei_Num and ((self.List_Num_Attached_Atom[I] == 3 and self.List_Num_Attached_Bond[I] == 3 ) or (self.List_Num_Attached_Atom[I] == 4 and self.List_Num_Attached_Bond[I] == 4 ))) :
                        Another_Subsidiary_Atom.append(I)
#                        print (I,"ANother_Subsidiary")
                    if (self.IAn[I] == self.C_Nuclei_Num and self.List_Num_Attached_Atom[I] == 3 and self.List_Num_Attached_Bond[I] == 4 ):
                        C_Atom.append(I)
#                        print (I,"C_Atom")
                if ( len(Central_Atom) == 1 and len(One_Subsidiary_Atom) == 1 and len(Another_Subsidiary_Atom) == 1 and len(C_Atom) == 2):
#                    print("Hello")
#                    if(self.BO[Central_Atom[0]][One_Subsidiary_Atom[0]] == 0 and self.BO[Central_Atom[0]][Another_Subsidiary_Atom[0]] == 0 and self.BO[Another_Subsidiary_Atom[0]][One_Subsidiary_Atom[0]] == 1 ):
                    if((self.BO[Central_Atom[0]][One_Subsidiary_Atom[0]] + self.BO[Central_Atom[0]][Another_Subsidiary_Atom[0]] + self.BO[Another_Subsidiary_Atom[0]][One_Subsidiary_Atom[0]]) == 1 ):
                        Num_Target_Group = Num_Target_Group + 1
                        for M in Ring:
                            self.Atom_Has_Been_Count.append(M)           #  adding the atom to the list of atom having been counted.
#                        Set_target_Group.append(Ring)
                Central_Atom = []
                One_Subsidiary_Atom = []
                Another_Subsidiary_Atom = []
                C_Atom = []
#                print(Num_Target_Group)
        return Num_Target_Group


    def Nine_Membered_Ring(self, Set_Phenyl_Group, Set_POX_Group):
        Num_Nine_Membered_Rings = 0
        for Ring_Phenyl in Set_Phenyl_Group:
            for Ring_POX in Set_POX_Group:
                identical_atom = list(set(Ring_Phenyl).intersection(set(Ring_POX)))
                if (len(identical_atom) == 2):
                    if ( self.BO[identical_atom[0]][identical_atom[1]] != 0):
                        Num_Nine_Membered_Rings = Num_Nine_Membered_Rings + 1
                        for M in Ring_Phenyl:
                            self.Atom_Has_Been_Count.append(M)           #  adding the atom to the list of atom having been counted.
                        for M in Ring_POX:
                            self.Atom_Has_Been_Count.append(M)           #  adding the atom to the list of atom having been counted.
        return Num_Nine_Membered_Rings

        
    
    def Five_Membered_Ring_With_Extra_Carbonyl(self):
        Num_C_Atom = 0
        Num_O_Atom = 0
        Num_carbonyl = 0
        Num_CO = 0
        Num_Five_Member_Ring_With_Extra_Carbonyl = 0
        for Ring in self.Ring_Set:
            Num_Ring = len(Ring)
            if (Num_Ring == 5):
                for I in Ring:
                    if (self.IAn[I] == self.C_Nuclei_Num):
                        Num_C_Atom = Num_C_Atom + 1
                    if (self.IAn[I] == self.O_Nuclei_Num):
                        Num_O_Atom = Num_O_Atom + 1
#                print ("O", Num_O_Atom, "C", Num_C_Atom)
                if (Num_O_Atom == 1 and Num_C_Atom == 4):
                    for I in Ring:
#                        print (self.IAn[I],self.List_Num_Attached_Atom[I], self.List_Num_Attached_Bond[I])
                        if (self.IAn[I] == self.C_Nuclei_Num and self.List_Num_Attached_Atom[I] == 3 and self.List_Num_Attached_Bond[I] == 4):
#                            print ("This is a message")
                            for J in self.List_Surrounding_Atom[I]:
                                if (self.IAn[J] == self.O_Nuclei_Num and self.BO[I][J] == 2 and (J not in Ring)):
                                    Num_carbonyl = Num_carbonyl + 1
                                    O_Outring_Mark = J # since an O atom is in the outside of ring, however, it is a part of the functional group, it must add to the Atom_Has_Been_Count.
                                if (self.IAn[J] == self.O_Nuclei_Num and self.BO[I][J] == 1 and (J in Ring)):
                                    Num_CO = Num_CO + 1
#                    print ("C_O", Num_CO, "C=O", Num_carbonyl)
                    if (Num_CO == 1 and Num_carbonyl == 1):
                        Num_Five_Member_Ring_With_Extra_Carbonyl = Num_Five_Member_Ring_With_Extra_Carbonyl + 1
                        self.Atom_Has_Been_Count.append(O_Outring_Mark)
                        for I in Ring:
                            self.Atom_Has_Been_Count.append(I)
                    Num_CO = 0
                    Num_carbonyl = 0
                Num_C_Atom = 0
                Num_O_Atom = 0
        return Num_Five_Member_Ring_With_Extra_Carbonyl

    def Five_Membered_Ring_With_Three_Coterminous_N_Atoms(self):
        Num_N_Atom = 0
        Num_C_Atom = 0
        Num_CC_Connected_Atom = 0
        Num_NN_Connected_Atom = 0
        Num_Functional_Group = 0
        for Ring in self.Ring_Set:
            if (len(Ring) == 5):
                for I in Ring:
                    if (self.IAn[I] == self.N_Nuclei_Num ):
                        Num_N_Atom = Num_N_Atom + 1
#                    if (self.IAn[I] == self.C_Nuclei_Num and self.List_Num_Attached_Atom[I] == 3 and self.List_Num_Attached_Bond == 4):
                    if (self.IAn[I] == self.C_Nuclei_Num):
                        Num_C_Atom = Num_C_Atom + 1
#                print ("Sign***C*", Num_C_Atom, "Sign****N***", Num_N_Atom, Ring)
                if (Num_C_Atom == 2 and Num_N_Atom == 3):
                    for I in Ring:
                        if (self.IAn[I] == self.C_Nuclei_Num):
                            for J in self.List_Surrounding_Atom[I]:
#                                print ("I, J, IAn[I],IAn[J],BO[I][J]", I, J, IAn[I], IAn[J], BO[I][J])
                                if (self.IAn[J] == self.C_Nuclei_Num and self.BO[I][J] == 2 and J in Ring):
                                    Num_CC_Connected_Atom = Num_CC_Connected_Atom + 1
#                        print("Num_CC_Double_Bond", Num_CC_Double_Bond)
                        if (self.IAn[I] == self.N_Nuclei_Num):
                            for J in self.List_Surrounding_Atom[I]:
                                if (self.IAn[J] == self.N_Nuclei_Num and self.BO[I][J] != 0 and J in Ring):
                                    Num_NN_Connected_Atom = Num_NN_Connected_Atom + 1
#                                if (self.IAn[J] == self.N_Nuclei_Num and self.BO[I][J] == 1):
#                                    Num_NN_Single_Bond = Num_NN_Single_Bond + 1
#                        print("Num_NN_Single_Bond, Num_NN_Double_Bond", Num_NN_Single_Bond, Num_NN_Double_Bond)
                    if (Num_NN_Connected_Atom/2 == 2 and Num_CC_Connected_Atom/2 == 1): 
                        Num_Functional_Group = Num_Functional_Group + 1
                        for I in Ring:
                            self.Atom_Has_Been_Count.append(I)
                Num_N_Atom = 0
                Num_C_Atom = 0
                Num_NN_Connected_Atom = 0
                Num_CC_Connected_Atom = 0
        return Num_Functional_Group



if __name__ == "__main__" :
    if (len(sys.argv) < 3) :
        print ("Usage: python Mutiple_Bond_Functional_Group.py BO Geom Category:1. cation 2. anion")
        sys.exit()

    Path1 = sys.argv[1]
    Path2 = sys.argv[2]
    cate  = sys.argv[3]
#    print("BO",os.path.getsize(Path1))
    if (os.path.getsize(Path1) == 0):
        print(Path2,"NO-BO")
        sys.exit()
#    IAn = np.loadtxt(Path2)
    BO = np.loadtxt(Path1)
    try:
        IAn = (np.loadtxt(Path2))[:,0]
    except IndexError:
        print (Path2, "Error,Halogen_ion")
        sys.exit()
    IAn = (np.loadtxt(Path2))[:,0]
#    IAn = (np.loadtxt(Path2))[:,0]
    ICoor = (np.loadtxt(Path2))[:,1:]
    NAtom = len(IAn)
    X = Mul_Ring(NAtom,BO,IAn,ICoor)
    _I_Cycle, _Final_Cycle = X.detect_cycle()



    #simple Functional Group 
#    __CH3  = X.Detecting_Simple_Connecting_Group(6,1,3)
#    __CH2  = X.Detecting_Simple_Connecting_Group(6,1,2)
#    __CH1  = X.Detecting_Simple_Connecting_Group(6,1,1)
#    __C_4  = X.Quarter_Atom(6,4)
#    __C_3  = X.Quarter_Atom(6,3)  #  carbocation
#    _Num_Halogen_Atom = X.Halogen_Counting()
#    __CF2_ = X.Detecting_Simple_Connecting_Group(6,9,2)
#    __CF3  = X.Detecting_Simple_Connecting_Group(6,9,3)
#    __COOC,_COOH,__COO_ = X.Detect_Ester()
#    __CO   = X.Double_Atom_Mutiple_Bond_Group(6,8,2,2)
#    __CHCH_= X.Double_Atom_Mutiple_Bond_Group(6,6,2,4,1,1)           #   Alkene
#    __CC_  = X.Double_Atom_Mutiple_Bond_Group(6,6,3,4)               #   Alkyne
#    __BH4  = X.Detecting_Simple_Connecting_Group(5,1,4)              #   _BH3
#    __B_4_ = X.Quarter_Atom(5,4)
#    __BF3  = X.Detecting_Simple_Connecting_Group(5,9,3)              #   _BF3
#    __O_   = X.Quarter_Atom(8,2)                                     #   Ether
#    __OH   = X.Detecting_Simple_Connecting_Group(8,1,1)              #   Hydroxyl
#    __NO2  = X.Detect_Nitro()                                        #   Nitro
#    __NH3  = X.Detecting_Simple_Connecting_Group(7,1,3)
#    __NH2  = X.Detecting_Simple_Connecting_Group(7,1,2)              #   amine
#    __NH1  = X.Detecting_Simple_Connecting_Group(7,1,1)              #   amine
#    __N_3  = X.Quarter_Atom(7,3)                                     #   Tertiary amine 
#    __CN   = X.Double_Atom_Mutiple_Bond_Group(6,7,3,3)               #   Cyano
#    _Set_Cyclohexane_Group, _Num_Cyclohexane_Group = X.Resemble_No_Hetero_Ring_Group(6,3)   # Cyclohexyl
#    _Num_O_Cyclopentyl_Group = X.Resemble_Single_Hetero_Ring_Group(5,3,8,2)                 # O_Cyclopentyl
#    _Set_Phenyl_Group, _Num_Phenyl_Group = X.Resemble_No_Hetero_Ring_Group(6,2)             # Phenyl
#    _Num_PO2_Group, _Set_PO2_Group = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,2,3,2,1,1)   # the number of PO2
#    __PO3_ = X.Detect_phosphate(1,2)                                 #   PO3
#    __PO4_ = X.Detect_phosphate(1,3)                                 #   PO4
#    __N_4_ = X.Quarter_Atom(7,4)                                     #   Quaternary ammonium
#    __N_2_ = X.Quarter_Atom(7,2)
#    __Al_4_ = X.Quarter_Atom(13,4)                                   #   Al+
#    __AlCl4_ = X.Detecting_Simple_Connecting_Group(13,17,4)            #   AlF4+
#    __P_4_ = X.Quarter_Atom(15,4)                                    #   P_4
#    __P_6_ = X.Quarter_Atom(15,6)                                    #   P_6
#    __S_3_ = X.Quarter_Atom(16,4)                                    #   Quaternary ammonium
#    __S_2_ = X.Quarter_Atom(16,2)                                    #   Quaternary ammonium
#    __Si_4_ = X.Quarter_Atom(14,4)                                    #   Quaternary ammonium
#    __SiH2_ = X.Detecting_Simple_Connecting_Group(14,1,2)            #   _SiH2_
#    __PF6_ = X.Detecting_Simple_Connecting_Group(15,9,6)             #   Hexafluorophosphate 
#    __SbF6_ = X.Detecting_Simple_Connecting_Group(51,9,6)            #   Fluoroantimonic 
#    __AsF6_ = X.Detecting_Simple_Connecting_Group(33,9,6)            #    
#    _Num_N_Phenys_Group = X.Resemble_Single_Hetero_Ring_Group(6,2,7,3)                      # N_Phenyl
#    _Num_Cyclopentyl_with_three_N = X.Five_Membered_Ring_With_Three_N_Atom()                # C2N3_
#    _Num_POX_Group_1, _Set_POX_Group_1 = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,3,3,2,2) # POX
#    _Num_POX_Group_2, _Set_POX_Group_2 = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,3,4,2,2) # POX
#    _Num_Nine_Membered_Rings_1 = X.Nine_Membered_Ring(_Set_Phenyl_Group, _Set_POX_Group_1)      # Nine membered ring
#    _Num_Nine_Membered_Rings_2 = X.Nine_Membered_Ring(_Set_Phenyl_Group, _Set_POX_Group_2)      # Nine membered ring
#    print("*******************************************")
#    _Num_POX_Anion, _Set_POX_Anion = X.Resemble_Double_Hetero_Ring_Group(5,7,2,2,7,2,2,2,2) # POX_Anion
#    print("*******************************************")
#    _Num_Nine_Membered_Rings_Anion = X.Nine_Membered_Ring(_Set_Phenyl_Group, _Set_POX_Anion)      # Nine membered ring for anion
#    print(_Num_POX_Anion, _Set_POX_Anion)
#    print(_Num_Nine_Membered_Rings_Anion)
#    _Num_PO0_Group, _Set_PO0_Group = X.Resemble_Double_Hetero_Ring_Group(6,7,4,4,8,2,2,3,3) # PO0  
#    __CN__ = X.Double_Atom_Mutiple_Bond_Group(6,7,2,4)               #   unamed 
#    _S_CN  = X.Detect_Thiocyanate()                                  #   Thiocyanate
#    _Num_N_Cyclopentyl_Group = X.Resemble_Single_Hetero_Ring_Group(5,3,7,4)                 # N_Cyclopentyl
#    _Num_PO1_Group, _Set_PO1_Group = X.Resemble_Double_Hetero_Ring_Group(6,7,3,3,7,2,2,2,2,0,1)            #   PO1  
#    _Num_PO3_Group, _Set_PO3_Group = X.Resemble_Double_Hetero_Ring_Group(5,16,3,3,7,2,3,2,2)               # PO3
#    _SO4_ = X.Detect_Sulfate(2,2)                                    # Sulfate _SO4_  
#    _SO3_ = X.Detect_Sulfate(2,1)                                    # Sulfate _SO3_  
#    _SO2_ = X.Detect_Sulfate(2,0)                                    # Sulfate _SO2_  
#    _Num_S_Cyclopentyl_Group = X.Resemble_Single_Hetero_Ring_Group(5,3,16,3)                # S_Cyclopentyl
    

    # Since some functional groups are calculated repeatedly, so there must be subtracted.
#    _Num_Halogen_Atom = _Num_Halogen_Atom-(2* __CF2_)-(3*__CF3)- 6*(__PF6_ + __AsF6_ + __SbF6_)
#    print(__CO,  __COO_, _COOH)
#    __CO = __CO - __COO_ - _COOH
#    __O_ = __O_ - __COOC - (2*__PO3_) 
#    __OH = __OH - _COOH
#    __CN  = __CN - _S_CN
#    _Num_Nine_Membered_Rings = _Num_Nine_Membered_Rings_1 + _Num_Nine_Membered_Rings_2
#    _Num_POX_Group = _Num_POX_Group_1 + _Num_POX_Group_2
#    _Num_Phenyl_Group = _Num_Phenyl_Group - _Num_Nine_Membered_Rings
#    _Num_POX_Group = _Num_POX_Group - _Num_Nine_Membered_Rings 
#    _Num_POX_Anion = _Num_POX_Anion - _Num_Nine_Membered_Rings_Anion 

#    print("helllllllllllllllllllllll")
#    print(X.BO[7][9],X.BO[7][8],__NO2)
#    print("helllllllllllllllllllllll")

#    Atoms_Havent_Been_Count = []
    
#    for I in X.List_Heavy_Atom:
#        if I not in list(set(X.Atom_Has_Been_Count)):
#            Atoms_Havent_Been_Count.append(I)
#    print(Atoms_Havent_Been_Count,"IIIIIIIIIII")
#    if (len(Atoms_Havent_Been_Count) != 0):
#        print(Path2,Atoms_Havent_Been_Count)
#        sys.exit()
#

#    print(__CH3, __CH2, __CH1, __C_4, _Num_Halogen_Atom, __CF2_, __CF3, __COOC, __COO_, _COOH, __CO, __CHCH_, __CC_, __O_, __OH, __NO2, __NH2, __NH1, __N_3, __CN, _Num_Cyclohexane_Group, _Num_O_Cyclopentyl_Group, _Num_Phenyl_Group, _Num_PO2_Group, __PO3_, __PO4_, __N_4_, __P_4_, __S_3_, _Num_N_Phenys_Group, _Num_Cyclopentyl_with_three_N, _Num_POX_Group, _Num_Nine_Membered_Rings, __NO2, _Num_PO0_Group, __CN__, _S_CN, _Num_N_Cyclopentyl_Group, _Num_PO1_Group, _Num_PO3_Group, _SO3_, _SO2_, _Num_S_Cyclopentyl_Group )
#    print(__CH3,"__CH3")
#    print(__CH2,"__CH2")
#    print(__CH1,"__CH1")
#    print(__C_4,"__C_4")
#    print(__C_3,"__C_3")
#    print(_Num_Halogen_Atom,"_Num_Halogen_Atom")
#    print(__CF2_,"__CF2_")
#    print(__CF3,"__CF3")
#    print(__COO_,"__COO_")
#    print(__COOC,"__COOC")
#    print(_COOH,"_COOH")
#    print(__CO,"__CO")
#    print(__CHCH_,"__CHCH_")
#    print( __CC_," __CC_")
#    print(__BH4,"__BH4")
#    print(__B_4_,"__B_4_")
#    print(__BF3,"__BF3")  
#    print(__O_,"__O_")
#    print(__OH,"__OH")
#    print(__NO2,"__NO2")
#    print(__NH2,"__NH2")
#    print(__NH1,"__NH1")
#    print(__N_3,"__N_3")
#    print(__CN,"__CN")
#    print(_Num_Cyclohexane_Group,"_Num_Cyclohexane_Group")
#    print(_Num_O_Cyclopentyl_Group,"_Num_O_Cyclopentyl_Group")
#    print(_Num_Phenyl_Group,"_Num_Phenyl_Group")
#    print(_Num_PO2_Group,"_Num_PO2_Group")
#    print(__PO3_,"__PO3_")
#    print(__PO4_,"__PO4_")
#    print(__N_4_,"__N_4_")
#    print(__N_2_,"__N_2_")
#    print( __P_4_," __P_4_")
#    print(__S_3_,"__S_3_")
#    print(__S_2_,"__S_2_")
#    print(__Si_4_,"__Si_4_")
#    print(__SiH2_,"__SiH2_")
#    print(__PF6_,"__PF6_")
#    print(__SbF6_,"__SbF6_")
#    print(__AsF6_,"__AsF6_")
#    print(_Num_N_Phenys_Group,"_Num_N_Phenys_Group")
#    print(_Num_Cyclopentyl_with_three_N,"_Num_Cyclopentyl_with_three_N")
#    print(_Num_POX_Group,"_Num_POX_Group")
#    print(_Num_Nine_Membered_Rings,"_Num_Nine_Membered_Rings")
#    print(__NO2,"__NO2")
#    print( _Num_PO0_Group," _Num_PO0_Group")
#    print(__CN__,"__CN__")
#    print(_S_CN,"_S_CN")
#    print(_Num_N_Cyclopentyl_Group,"_Num_N_Cyclopentyl_Group")
#    print(_Num_PO1_Group,"_Num_PO1_Group")
#    print(_Num_PO3_Group,"_Num_PO3_Group")
#    print(_SO4_,"_SO4_")
#    print(_SO3_,"_SO3_")
#    print(_SO2_,"_SO2_")
#    print(_Num_S_Cyclopentyl_Group,"_Num_S_Cyclopentyl_Group")



# for cation 
#    print(cate)
    if (int(cate) == 1):
        __CH3  = X.Detecting_Simple_Connecting_Group(6,1,3)
        __CH2  = X.Detecting_Simple_Connecting_Group(6,1,2)
        __CH1  = X.Detecting_Simple_Connecting_Group(6,1,1)
        __C_4  = X.Quarter_Atom(6,4)
        __C_3  = X.Quarter_Atom(6,3)
        _Num_Halogen_Atom = X.Halogen_Counting()
        __CF2_ = X.Detecting_Simple_Connecting_Group(6,9,2)
        __CF3  = X.Detecting_Simple_Connecting_Group(6,9,3)
        __COOC,_COOH,__COO_ = X.Detect_Ester()
        __CHCH_= X.Double_Atom_Mutiple_Bond_Group(6,6,2,4,1,1)           #   Alkene
        __CC_  = X.Double_Atom_Mutiple_Bond_Group(6,6,3,4)               #   Alkyne
        __O_   = X.Quarter_Atom(8,2)                                     #   Ether
        __OH   = X.Detecting_Simple_Connecting_Group(8,1,1)              #   Hydroxyl
        __CO   = X.Double_Atom_Mutiple_Bond_Group(6,8,2,2)
        __NO2  = X.Detect_Nitro()                                        #   Nitro
        __NH3  = X.Detecting_Simple_Connecting_Group(7,1,3)
        __NH2  = X.Detecting_Simple_Connecting_Group(7,1,2)              #   amine
        __NH1  = X.Detecting_Simple_Connecting_Group(7,1,1)              #   amine
        __N_3  = X.Quarter_Atom(7,3)                                     #   Tertiary amine 
        __N_4_ = X.Quarter_Atom(7,4)                                     #   Quaternary ammonium
        __CN   = X.Double_Atom_Mutiple_Bond_Group(6,7,3,3)               #   Cyano
        _Set_Cyclohexane_Group, _Num_Cyclohexane_Group = X.Resemble_No_Hetero_Ring_Group(6,3)   # Cyclohexane
        _Set_Cyclopentane_Group, _Num_Cyclopentane_Group = X.Resemble_No_Hetero_Ring_Group(5,3)   # Cyclopentane
        _Num_O_Cyclopentyl_Group = X.Resemble_Single_Hetero_Ring_Group(5,3,8,2)                 # O_Cyclopentyl
        _Set_Phenyl_Group, _Num_Phenyl_Group = X.Resemble_No_Hetero_Ring_Group(6,2)             # Phenyl
        _Num_PO2_Group_1, _Set_PO2_Group_1 = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,2,3,2,1,1)   # the number of PO2
        _Num_PO2_Group_2, _Set_PO2_Group_2 = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,3,4,2,1,1)   # the number of PO2
        _Num_PO2_Group = _Num_PO2_Group_2 + _Num_PO2_Group_1
        __PO3_ = X.Detect_phosphate(1,2)                                 #   PO3
        __PO4_ = X.Detect_phosphate(1,3)                                 #   PO4
        __P_4_ = X.Quarter_Atom(15,4)                                    #   P_4
        __Si_4_ = X.Quarter_Atom(14,4)                                    #   Quaternary ammonium
        __SiH2_ = X.Detecting_Simple_Connecting_Group(14,1,2)            #   _SiH2_
        _Num_N_Phenys_Group_1 = X.Resemble_Single_Hetero_Ring_Group(6,2,7,3)                      # N_Phenyl
        _Num_N_Phenys_Group_2 = X.Resemble_Single_Hetero_Ring_Group(6,2,7,2)                      # N_Phenyl, in this case, N conncet wth 2 atoms
        _Num_Cyclopentyl_with_three_N = X.Five_Membered_Ring_With_Three_N_Atom()                # C2N3_
        _Num_POX_Group_1, _Set_POX_Group_1 = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,3,3,2,2) # POX
        _Num_POX_Group_2, _Set_POX_Group_2 = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,3,4,2,2) # POX
 #       print("XXXXXXXXXXXXXXXXXXX")
        _Num_POX_Group_3, _Set_POX_Group_3 = X.Resemble_Double_Hetero_Ring_Group(5,7,4,4,7,2,3,2,2) # POX
        _Num_POX_Group_4, _Set_POX_Group_4 = X.Resemble_Double_Hetero_Ring_Group(5,7,4,4,7,3,4,2,2) # POX
 #       print("XXXXXXXXXXXXXXXXXXX")
        _Num_Nine_Membered_Rings_1 = X.Nine_Membered_Ring(_Set_Phenyl_Group, _Set_POX_Group_1)      # Nine membered ring
        _Num_Nine_Membered_Rings_2 = X.Nine_Membered_Ring(_Set_Phenyl_Group, _Set_POX_Group_2)      # Nine membered ring
        _Num_PO0_Group, _Set_PO0_Group = X.Resemble_Double_Hetero_Ring_Group(6,7,4,4,8,2,2,3,3) # PO0
        __CN__1 = X.Double_Atom_Mutiple_Bond_Group(6,7,2,4)               #   __CN__
        __CN__2 = X.Double_Atom_Mutiple_Bond_Group(6,7,2,3)               #   __CN__ in this case, N only conncet with 2 aotms in which one is single bond another is double bond.
        _Num_N_Cyclopentyl_Group = X.Resemble_Single_Hetero_Ring_Group(5,3,7,4)                 # N_Cyclopentyl
        _Num_N_Cyclhexane_Group = X.Resemble_Single_Hetero_Ring_Group(6,3,7,4)                 # N_Cyclhexane
        _Num_PO1_Group, _Set_PO1_Group = X.Resemble_Double_Hetero_Ring_Group(6,7,2,2,7,3,3,2,2)            #   PO1
        _Num_PO3_Group, _Set_PO3_Group = X.Resemble_Double_Hetero_Ring_Group(5,16,3,3,7,2,3,2,2)               # PO3
        _Num_POA_Group, _Set_POA_Group = X.Resemble_Double_Hetero_Ring_Group(5,7,3,4,7,3,3,2,2)               # POA
        _Num_S_Cyclopentyl_Group = X.Resemble_Single_Hetero_Ring_Group(5,3,16,3)                # S_Cyclopentyl
        _Num_S_Cyclohexane_Group = X.Resemble_Single_Hetero_Ring_Group(6,3,16,3)                # S_Cyclohexane
        _Num_Five_Member_Ring_With_Extra_Carbonyl = X.Five_Membered_Ring_With_Extra_Carbonyl()
        _Num_Five_Membered_Ring_With_Three_Coterminous_N_Atoms = X.Five_Membered_Ring_With_Three_Coterminous_N_Atoms()
        __S_3_ = X.Quarter_Atom(16,3)                                    #   Quaternary ammonium
        __S_2_ = X.Quarter_Atom(16,2)                                    #   -S-
                                                                         #      |    
        __S_6_ = X.Quarter_Atom(16,6)                                    #     =S=
                                                                         #      |        
        _SO2_ = X.Detect_Sulfate(2,0)
        _SO3_ = X.Detect_Sulfate(2,1)
        _SO4_ = X.Detect_Sulfate(2,2)                                    # Sulfate _SO4_
        _SF5_ = X.Detecting_Simple_Connecting_Group(16,9,5)
        __Si_4_ = X.Quarter_Atom(14,4)                                    #   Quaternary ammonium
        __SiH2_ = X.Detecting_Simple_Connecting_Group(14,1,2)            #   _SiH2_
        __SiH3  = X.Detecting_Simple_Connecting_Group(14,1,3)            #   _SiH3
#        print(X.BO[1][6],X.List_Num_Attached_Atom[6],X.List_Num_Attached_Bond[6])
        # Since some functional groups are calculated repeatedly, so there must be subtracted.
        _Num_Halogen_Atom = _Num_Halogen_Atom-(2* __CF2_)-(3*__CF3)
        __CO = __CO  - _COOH - __COOC
        __O_ = __O_ - __COOC - (1*__PO3_) - (2* __PO4_) 
        __OH = __OH - _COOH
        __CN__ = __CN__2 + __CN__1
#        __CN  = __CN - _S_C
        _Num_N_Phenys_Group = _Num_N_Phenys_Group_1 + _Num_N_Phenys_Group_2
        _Num_Nine_Membered_Rings = _Num_Nine_Membered_Rings_1 + _Num_Nine_Membered_Rings_2
        _Num_POX_Group = _Num_POX_Group_1 + _Num_POX_Group_2+ _Num_POX_Group_3 + _Num_POX_Group_4
        _Num_Phenyl_Group = _Num_Phenyl_Group - _Num_Nine_Membered_Rings
        _Num_POX_Group = _Num_POX_Group - _Num_Nine_Membered_Rings 
#        print(X.BO[4])
#        print(X.BO[39])
#        print(X.BO[2])
#        print(X.BO[3])
#        print(X.BO[4])

#        print (_Num_Five_Membered_Ring_With_Three_Coterminous_N_Atoms)
#        print(X.BO[0][9],X.BO[5][9])
#        print(X.BO[9])
#        print(X.BO[7])
        Atoms_Havent_Been_Count = []
#         print(X.BO[0][1])
#        print (X.Atom_Has_Been_Count)
#        print(X.BO[13])
#        print(X.BO[1])
#        print(X.BO[2])
#        print(X.BO[3])
#        print(X.BO[4])
#        print(X.BO[64])
#        print(X.BO[32])
#        print(X.List_Circle_Atom)
#        print(X.Ring_Set)
        for I in X.List_Heavy_Atom:
            if I not in list(set(X.Atom_Has_Been_Count)):
                Atoms_Havent_Been_Count.append(I)
#                print(Atoms_Havent_Been_Count,"IIIIIIIIIII")
#                print(BO[49])
        if (len(Atoms_Havent_Been_Count) != 0):
            #        if (len(Atoms_Havent_Been_Count) == 6)
#            print(Path2,Atoms_Havent_Been_Count)
            sys.exit()
#        print(__NH1)
        print(Path2, _Num_Cyclopentane_Group)
#        print(Path2, __CH3, __CH2, __CH1, __CO, __COO_ , _COOH,  __COOC, __NH1,__NH2, __NH3,  __N_4_, __O_, _Num_Halogen_Atom, __OH, __C_4, __CHCH_, __CC_, __P_4_,__CF3, __CF2_, __NO2, __S_2_, __S_3_, __N_3, __Si_4_, __S_6_, __CN, __CN__, _Num_Phenyl_Group, _Num_Cyclopentyl_with_three_N, _Num_POX_Group, _Num_N_Phenys_Group, _Num_Cyclohexane_Group, _Num_N_Cyclhexane_Group, _Num_Five_Membered_Ring_With_Three_Coterminous_N_Atoms, _Num_S_Cyclohexane_Group, _Num_Five_Member_Ring_With_Extra_Carbonyl, _Num_Nine_Membered_Rings, _Num_PO2_Group, _Num_N_Cyclopentyl_Group, _Num_S_Cyclopentyl_Group, _Num_PO0_Group , _Num_O_Cyclopentyl_Group, __PO3_,  _Num_PO1_Group, _Num_PO3_Group, _SO3_ )


# for anion
    if (int(cate) == 2):
#
        '''
        __CH3  = X.Detecting_Simple_Connecting_Group(6,1,3)
        __CH2  = X.Detecting_Simple_Connecting_Group(6,1,2)
        __CO   = X.Double_Atom_Mutiple_Bond_Group(6,8,2,2)
        __NO2  = X.Detect_Nitro()                                        #   Nitro
        __B_4_ = X.Quarter_Atom(5,4)
        __NO3  = X.Detect_Nitrate()
        __COOC,_COOH,__COO_ = X.Detect_Ester()
        _Num_Halogen_Atom = X.Halogen_Counting()
        _SO2_ = X.Detect_Sulfate(2,0)
        _SO4_ = X.Detect_Sulfate(2,2)                                    # Sulfate _SO4_
        __CH1  = X.Detecting_Simple_Connecting_Group(6,1,1)
        __NH2  = X.Detecting_Simple_Connecting_Group(7,1,2)
        __NH1  = X.Detecting_Simple_Connecting_Group(7,1,1)
        __C_4 = X.Quarter_Atom(6,4)
        __CHCH_= X.Double_Atom_Mutiple_Bond_Group(6,6,2,4,1,1)           #   Alkene
        __CN   = X.Double_Atom_Mutiple_Bond_Group(6,7,3,3)               #   Cyano
        __N_3  = X.Quarter_Atom(7,3)                                     #   Tertiary amine
        __CF3  = X.Detecting_Simple_Connecting_Group(6,9,3)
        __CF2_ = X.Detecting_Simple_Connecting_Group(6,9,2)
        __BF3_ = X.Detecting_Simple_Connecting_Group(5,9,3)
        __BF4_ = X.Detecting_Simple_Connecting_Group(5,9,4)
        __BH4_ = X.Detecting_Simple_Connecting_Group(5,1,4)
        __PO4_ = X.Detect_phosphate(1,3)                                 #   PO4
        __N_2_ = X.Quarter_Atom(7,2)
        __S_2_ = X.Quarter_Atom(16,2)                                    #   Quaternary ammonium
        __SH  = X.Detecting_Simple_Connecting_Group(16,1,1)
        __SiH2_ = X.Detecting_Simple_Connecting_Group(14,1,2)            #   _SiH2_
        __O_   = X.Quarter_Atom(8,2)                                     #   Ether
        __OH   = X.Detecting_Simple_Connecting_Group(8,1,1)              #   Hydroxyl
        _Num_POX_Group_1, _Set_POX_Group_1 = X.Resemble_Double_Hetero_Ring_Group(5,7,2,2,7,2,2,2,2) # POX
        _Num_POX_Group_2, _Set_POX_Group_2 = X.Resemble_Double_Hetero_Ring_Group(5,7,2,2,7,2,3,2,2) # POX
        _Num_POX_Group_3, _Set_POX_Group_3 = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,2,3,2,2) # POX
        _Num_POX_Group_4, _Set_POX_Group_4 = X.Resemble_Double_Hetero_Ring_Group(5,7,4,4,7,3,4,2,2) # POX
        __SbF6_ = X.Detecting_Simple_Connecting_Group(51,9,6)            #   Fluoroantimonic
        __Al_4_ = X.Quarter_Atom(13,4)                                   #   Al+
        __Fe_4_ = X.Quarter_Atom(56,4)                                   #   Fe+
#        __P_6_ = X.Quarter_Atom(15,6)                                   #   P6+
        _S_CN  = X.Detect_Thiocyanate()                                  #   Thiocyanate
#       N=N=N
        _N_C_N  = X.N_C_N()                                              #   Cyanamide
        _Set_Cyclohexane_Group, _Num_Cyclohexane_Group = X.Resemble_No_Hetero_Ring_Group(6,3)   # Cyclohexane
        _Num_PO0_Group, _Set_PO0_Group = X.Resemble_Double_Hetero_Ring_Group(6,7,4,4,8,2,2,3,3) # PO0
#       phencyl with O atom 
        _Num_Five_Membered_Ring_With_Three_Coterminous_N_Atoms = X.Five_Membered_Ring_With_Three_Coterminous_N_Atoms()
#       Halogen ion
        _SO3_ = X.Detect_Sulfate(2,1)
        _Num_PO2_Group_1, _Set_PO2_Group_1 = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,2,3,2,1,1)   # the number of PO2
        _Num_PO2_Group_2, _Set_PO2_Group_2 = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,3,4,2,1,1)   # the number of PO2
        _Num_PO2_Group = _Num_PO2_Group_2 + _Num_PO2_Group_1
        __PO2_ = X.Detect_phosphate(1,1)                                 #   PO2
        __PO3_ = X.Detect_phosphate(1,2)                                 #   PO3
        __C_3  = X.Quarter_Atom(6,3)  #  carbocation
        _Set_Phenyl_Group, _Num_Phenyl_Group = X.Resemble_No_Hetero_Ring_Group(6,2) # Phenyl
        _Num_Cyclopentyl_with_three_N = X.Five_Membered_Ring_With_Three_N_Atom()                # C2N3_
        __AsF6_ = X.Detecting_Simple_Connecting_Group(33,9,6)            #
        _Num_Nine_Membered_Rings_1 = X.Nine_Membered_Ring(_Set_Phenyl_Group, _Set_POX_Group_1)      # Nine membered ring
        _Num_Nine_Membered_Rings_2 = X.Nine_Membered_Ring(_Set_Phenyl_Group, _Set_POX_Group_2)      # Nine membered ring
        _Num_Nine_Membered_Rings_3 = X.Nine_Membered_Ring(_Set_Phenyl_Group, _Set_POX_Group_3)      # Nine membered ring
        _Num_Nine_Membered_Rings = _Num_Nine_Membered_Rings_1 + _Num_Nine_Membered_Rings_2 + _Num_Nine_Membered_Rings_3
        __PF6_ = X.Detecting_Simple_Connecting_Group(15,9,6)             #   Hexafluorophosphate
        _Num_N_Cyclopentyl_Group_1 = X.Resemble_Single_Hetero_Ring_Group(5,3,7,4)                 # N_Cyclopentyl
        _Num_N_Cyclopentyl_Group_2 = X.Resemble_Single_Hetero_Ring_Group(5,3,7,3)                 # N_Cyclopentyl
        _Num_N_Cyclopentyl_Group = _Num_N_Cyclopentyl_Group_1 + _Num_N_Cyclopentyl_Group_2
        _Num_S_Cyclopentadiene_Group = X.Resemble_Single_Hetero_Ring_Group(5,2,16,2)                # S_Cyclopentadiene
        _Num_O_Cyclopentadiene_Group = X.Resemble_Single_Hetero_Ring_Group(5,2,8,2)                # O_Cyclopentadiene
        _Num_PO2_Group_1, _Set_PO2_Group_1 = X.Resemble_Double_Hetero_Ring_Group(5,7,2,2,7,2,3,2,1,1)   # the number of PO2
        _Num_PO2_Group_2, _Set_PO2_Group_2 = X.Resemble_Double_Hetero_Ring_Group(5,7,2,2,7,2,2,2,1,1)   # the number of PO2
        _Num_PO2_Group = _Num_PO2_Group_2 + _Num_PO2_Group_1

        __CO = __CO - __COO_ - _COOH - __COOC
        __O_ = __O_ - __COOC - (1*__PO3_) - (2* __PO4_)  - (1* __PO2_) - (1* _SO3_) - (2* _SO4_ ) - (1* __COOC) 
        __OH = __OH - _COOH
        __CN  = __CN - _S_CN

        _Num_Halogen_Atom = _Num_Halogen_Atom - (3*__CF3) - (2*__CF2_) - 6* (__PF6_ + __SbF6_ + __AsF6_) - 4 *( __BF4_) -3 *(__BF3_)
        _Num_Phenyl_Group = _Num_Phenyl_Group - _Num_Nine_Membered_Rings
        _Num_POX_Group = _Num_POX_Group_1 + _Num_POX_Group_2+ _Num_POX_Group_3 + _Num_POX_Group_4 - _Num_Nine_Membered_Rings
        '''


        __CH3  = X.Detecting_Simple_Connecting_Group(6,1,3)
        __CH2  = X.Detecting_Simple_Connecting_Group(6,1,2)
        __CH1  = X.Detecting_Simple_Connecting_Group(6,1,1)
        __C_4 = X.Quarter_Atom(6,4)
        __C_3  = X.Quarter_Atom(6,3)  #  carbocation
        __CHCH_= X.Double_Atom_Mutiple_Bond_Group(6,6,2,4,1,1)           #   Alkene
        __CO   = X.Double_Atom_Mutiple_Bond_Group(6,8,2,2)
        __CN   = X.Double_Atom_Mutiple_Bond_Group(6,7,3,3)               #   Cyano
        _S_CN  = X.Detect_Thiocyanate()                                  #   Thiocyanate
        _N_C_N  = X.N_C_N()                                              #   Cyanamide
        __CF3  = X.Detecting_Simple_Connecting_Group(6,9,3)
        __CF2_ = X.Detecting_Simple_Connecting_Group(6,9,2)

        __NO3  = X.Detect_Nitrate()
        __NO2  = X.Detect_Nitro()                                        #   Nitro
        __NH2  = X.Detecting_Simple_Connecting_Group(7,1,2)
        __NH1  = X.Detecting_Simple_Connecting_Group(7,1,1)
        __N_3  = X.Quarter_Atom(7,3)                                     #   Tertiary amine
        __N_2_ = X.Quarter_Atom(7,2)
        _Num_Halogen_Atom = X.Halogen_Counting()

        _SO2_ = X.Detect_Sulfate(2,0)
        _SO3_ = X.Detect_Sulfate(2,1)
        _SO3H_ = X.Detect_Sulfate(2,1,1)
        _SO3_Total = _SO3_ + _SO3H_ 
        _SO4_ = X.Detect_Sulfate(2,2)                                    # Sulfate _SO4_
        _SO4H_ = X.Detect_Sulfate(2,2,1)                                    # Sulfate _SO4H_
        _SO4_Total = _SO4H_ + _SO4_ 
        __S_2_ = X.Quarter_Atom(16,2)                                    #   Quaternary ammonium
        __SH  = X.Detecting_Simple_Connecting_Group(16,1,1)
        __SiH2_ = X.Detecting_Simple_Connecting_Group(14,1,2)            #   _SiH2_
        __Si_4_ = X.Quarter_Atom(14,4)                                    #   Quaternary ammonium
        __O_   = X.Quarter_Atom(8,2)                                     #   Ether
        __OH   = X.Detecting_Simple_Connecting_Group(8,1,1)              #   Hydroxyl
        __COOC,_COOH,__COO_ = X.Detect_Ester()
        __PO2_ = X.Detect_phosphate(1,1)                                 #   PO2
        __PO3_ = X.Detect_phosphate(1,2)                                 #   PO3
        __PO3H_ = X.Detect_phosphate(1,2,1)                                 #   PO3H
        __PO3_Total = __PO3_ + __PO3H_
        __PO4H2_ = X.Detect_phosphate(1,3,2)                                 #   PO4
        __PO4_ = X.Detect_phosphate(1,3)                                 #   PO4H2
        __PO4_Total = __PO4_ + __PO4H2_
        __B_4_ = X.Quarter_Atom(5,4)
        __BF3_ = X.Detecting_Simple_Connecting_Group(5,9,3)
        __BF4_ = X.Detecting_Simple_Connecting_Group(5,9,4)
        __BH4_ = X.Detecting_Simple_Connecting_Group(5,1,4)
        __Al_4_ = X.Quarter_Atom(13,4)                                   #   Al+
        __Fe_4_ = X.Quarter_Atom(56,4)                                   #   Fe+
        __SbF6_ = X.Detecting_Simple_Connecting_Group(51,9,6)            #   Fluoroantimonic
        __AsF6_ = X.Detecting_Simple_Connecting_Group(33,9,6)            #
        __PF6_ = X.Detecting_Simple_Connecting_Group(15,9,6)             #   Hexafluorophosphate
        _Num_POX_Group_1, _Set_POX_Group_1 = X.Resemble_Double_Hetero_Ring_Group(5,7,2,2,7,2,2,2,2) # POX
        _Num_POX_Group_2, _Set_POX_Group_2 = X.Resemble_Double_Hetero_Ring_Group(5,7,2,2,7,2,3,2,2) # POX
        _Num_POX_Group_3, _Set_POX_Group_3 = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,2,3,2,2) # POX
        _Num_POX_Group_4, _Set_POX_Group_4 = X.Resemble_Double_Hetero_Ring_Group(5,7,4,4,7,3,4,2,2) # POX

        _Set_Cyclohexane_Group, _Num_Cyclohexane_Group = X.Resemble_No_Hetero_Ring_Group(6,3)   # Cyclohexane
        _Num_PO0_Group, _Set_PO0_Group = X.Resemble_Double_Hetero_Ring_Group(6,7,4,4,8,2,2,3,3) # PO0
        _Num_Five_Membered_Ring_With_Three_Coterminous_N_Atoms = X.Five_Membered_Ring_With_Three_Coterminous_N_Atoms()
        _Num_PO2_Group_1, _Set_PO2_Group_1 = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,2,3,2,1,1)   # the number of PO2
        _Num_PO2_Group_2, _Set_PO2_Group_2 = X.Resemble_Double_Hetero_Ring_Group(5,7,3,3,7,3,4,2,1,1)   # the number of PO2
        _Num_PO2_Group_3, _Set_PO2_Group_3 = X.Resemble_Double_Hetero_Ring_Group(5,7,2,3,7,2,2,2,1,1)   # the number of PO2
        _Num_PO2_Group = _Num_PO2_Group_2 + _Num_PO2_Group_1 + _Num_PO2_Group_3
        _Set_Phenyl_Group, _Num_Phenyl_Group = X.Resemble_No_Hetero_Ring_Group(6,2) # Phenyl
        _Num_Cyclopentyl_with_three_N = X.Five_Membered_Ring_With_Three_N_Atom()                # C2N3_

        _Num_Nine_Membered_Rings_1 = X.Nine_Membered_Ring(_Set_Phenyl_Group, _Set_POX_Group_1)      # Nine membered ring
        _Num_Nine_Membered_Rings_2 = X.Nine_Membered_Ring(_Set_Phenyl_Group, _Set_POX_Group_2)      # Nine membered ring
        _Num_Nine_Membered_Rings_3 = X.Nine_Membered_Ring(_Set_Phenyl_Group, _Set_POX_Group_3)      # Nine membered ring
        _Num_Nine_Membered_Rings = _Num_Nine_Membered_Rings_1 + _Num_Nine_Membered_Rings_2 + _Num_Nine_Membered_Rings_3

        _Num_S_Cyclopentadiene_Group = X.Resemble_Single_Hetero_Ring_Group(5,2,16,2)                # S_Cyclopentadiene
        _Num_O_Cyclopentadiene_Group = X.Resemble_Single_Hetero_Ring_Group(5,2,8,2)

        __CO = __CO - __COO_ - _COOH - __COOC
        __O_ = __O_ - (1*__PO3_) - (2* __PO4_) - (1* _SO4_ ) - (1* __COOC)
        __OH = __OH - _COOH -  2 * __PO4H2_ - 1 * __PO3H_ - 1 * _SO4H_  - _SO3H_ 
        __CN  = __CN - _S_CN

        _Num_Phenyl_Group = _Num_Phenyl_Group - _Num_Nine_Membered_Rings
        _Num_POX_Group = _Num_POX_Group_1 + _Num_POX_Group_2+ _Num_POX_Group_3 + _Num_POX_Group_4 - _Num_Nine_Membered_Rings
        _Num_Halogen_Atom = _Num_Halogen_Atom - (3 * __CF3) - (2 * __CF2_) - (3 * __BF3_) - (4 * __BF4_) - 6 * (__SbF6_ + __AsF6_ + __PF6_)


        Atoms_Havent_Been_Count = []
#        print(X.List_Surrounding_Atom[2][1:])
#        print(X.BO[0])
#        print( __COOC, _COOH, __COO_)
        for I in X.List_Heavy_Atom:
            if I not in list(set(X.Atom_Has_Been_Count)):
                Atoms_Havent_Been_Count.append(I)
#        print(Atoms_Havent_Been_Count,"IIIIIIIIIII")
        if (len(Atoms_Havent_Been_Count) != 0):
#            print(Path2,Atoms_Havent_Been_Count)
#            print(BO[8])
            sys.exit()
#        print ("1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 ")
#        print( __COOC, _COOH, __COO_)
        print(Path2, __CH3, __CH2, __CH1, __C_4, __C_3, __CHCH_, __CO, __CN ,  _S_CN, _N_C_N, __CF3, __CF2_, __NO3, __NO2, __NH2, __NH1, __N_3, __N_2_, _Num_Halogen_Atom, _SO2_, _SO3_Total, _SO4_Total, __S_2_, __SH, __SiH2_, __Si_4_, __O_, __OH, _COOH, __COOC, __COO_, __PO2_, __PO3_Total, __PO4_Total, __B_4_, __BF3_, __BF4_, __BH4_, __Al_4_, __Fe_4_, __SbF6_, __AsF6_, __PF6_, _Num_POX_Group, _Num_Cyclohexane_Group, _Num_PO0_Group, _Num_Five_Membered_Ring_With_Three_Coterminous_N_Atoms, _Num_PO2_Group, _Num_Phenyl_Group, _Num_Cyclopentyl_with_three_N,  _Num_Nine_Membered_Rings, _Num_S_Cyclopentadiene_Group,  _Num_O_Cyclopentadiene_Group)


