import numpy as np 
import sys
sys.path.append('../')
import os
#from Simple_INFO import Single_Bond

#try:
 #   from .Simple_INFO import Single_Bond
#except Exception: #ImportError
  #  from Simple_INFO import Single_Bond
from Cycle_Detecting import Ring_Detecting


#class Group_Single_Bond(Single_Bond,Ring_Detecting):
class Group_Single_Bond(Ring_Detecting):
    def __init__(self,NAtom,BO,IAn,ICoor): # Nomal_Bond is the number of bond that the initial atom aonnected with other, for instance, C is 4,O is 2
#        Single_Bond.__init__(self,NAtom,BO,IAn,ICoor)
        Ring_Detecting.__init__(self,NAtom,BO,IAn,ICoor)
#        self.Initial_Atom                                       = Initial_Atom
#        self.List_Num_Attached_Atom,self.List_Num_Attached_Bond = Single_Bond.Attach(self)
#        self.List_Surrounding_Atom                              = Single_Bond.Surrounding_Atom(self)
#        self.List_Hybrid                                        = Single_Bond.Hybrid(self,self.List_Num_Attached_Bond,self.List_Num_Attached_Atom,self.List_Surrounding_Atom)
#        self.List_Num_of_terminal_atom_Surrounding              = Single_Bond.Num_of_terminal_atom_Surrounding(self)
        self.List_Num_Attached_Atom,self.List_Num_Attached_Bond = Ring_Detecting.Attach(self)
        self.List_Surrounding_Atom                              = Ring_Detecting.Surrounding_Atom(self)
        self.List_Hybrid                                        = Ring_Detecting.Hybrid(self,self.List_Num_Attached_Bond,self.List_Num_Attached_Atom,self.List_Surrounding_Atom)
        self.List_Num_of_terminal_atom_Surrounding              = Ring_Detecting.Num_of_terminal_atom_Surrounding(self)
        self.Ring_Set, self.List_Circle_Atom                    = Ring_Detecting.detect_cycle(self)
#       **********************


    
    def Terminal_Group(self): # this method is aim at detecting the Terminal group like -CH3,-NH2,-CF2H,-OH
        M = 0 
        Icount = 0
        List_Num_Halogen_and_H = []
        List_Terminal_group_Central_Atom = []
        for I in self.List_Surrounding_Atom:    
            for J in I[1:] :
                if (self.BO[I[0],J] == 1 and self.IAn[J] in self.Terminal_Atom):  # if the attached atom is H or halogen
                    Icount = Icount + 1
            List_Num_Halogen_and_H.append(Icount)                            
            Icount = 0
        while M < self.NAtom:
            if (self.List_Hybrid[M] == 3 ):                                       # if the central atom is SP3 hybirdration, the sp and sp2 hybirdration are not taken into consideration
                if (self.IAn[M] == self.C_Nuclei_Num and List_Num_Halogen_and_H[M] ==3  ):  # if the central atom is C, and the number of connected atom which is Halogen or H is 3, which means that the C is a terminal C
                    List_Terminal_group_Central_Atom.append(1)
                    M = M + 1
                    continue
     
                    M = M + 1
                    continue
                if (self.IAn[M] == self.O_Nuclei_Num and List_Num_Halogen_and_H[M] ==1  ):  # if the central atom is O, and the number of connected atom which is Halogen or H is 1, which means that the O is a terminal O
                    List_Terminal_group_Central_Atom.append(1)
                    M = M + 1
                    continue
            List_Terminal_group_Central_Atom.append(0)
            M = M + 1    
        return List_Terminal_group_Central_Atom   
        

    
    def Detecting_Simple_Connecting_Group(self,Central_Atom_Type,Attaching_Atom_Type,Num_Attaching_Atom): # eg. -CH3 (self,6,1,3)  -CF2- (self,6, 9, 2) -NH2 (7,1,2), can detect -CH3,-CH2-,-CH-,-NH2,-NH-,-OH,-CF3,-CF2-.
        I = 0
        J = 0
        ICount = 0
        Group_num = 0
        while I < self.NAtom :
#            if (self.IAn[I] == 5):
#                print(I,self.IAn[I], self.List_Circle_Atom[I],"BBBBBBBBBBBBBBBBBBBBBB")
            if (self.IAn[I] == Central_Atom_Type and 2 not in self.BO[I] and 3 not in self.BO[I] and self.List_Circle_Atom[I] == 0):
                while J < self.NAtom:
                    if (self.BO[I][J] == 1 and self.IAn[J] == Attaching_Atom_Type ):  # if J Atom connect with C with single bond and J atom nuclei number is attached_Atom_Type
                        ICount = ICount + 1
                    J =J + 1
#                print(I,self.IAn[I],Attaching_Atom_Type,ICount,"helloooooo")
                if (ICount == Num_Attaching_Atom ):   # the number of one kind of atom connected with C in single bond, for instance -CH3 AttachedAtom_Number is 3 
                    Group_num = Group_num + 1
                    self.Atom_Has_Been_Count.append(I)
                ICount = 0
            I = I + 1
            J = 0  # when the initial atom changes, J has to reset to 0
        return Group_num



    def Halogen_Counting(self):       # the method is aim at counting the number of Halogen Atom(F,Cl,Br) in the whole molecular or ion
        I = 0
        Number_Halogen_Atom = 0
        while I < self.NAtom:
            if (self.IAn[I] in self.Halogen_Atom):
                Number_Halogen_Atom = Number_Halogen_Atom + 1 
                self.Atom_Has_Been_Count.append(I)           #  adding the atom to the list of atom having been counted.  
            I = I + 1
        return Number_Halogen_Atom 
    


    def Quarter_Atom(self,Central_Atom_Type,Num_Attaching_C_Atom):         #                                                               |     |       |    |            |            |
  # This method is designed for detecting the Quarter atom in molecular, like C N S P ,all connecting bond are single bond, This method can detect -C- , -P(+)-, -N-, -N(+)-, -S-, -S(-)-, -O-, -B(+)-
  #                                                                                                                                                 |     |            |                         |
        I = 0 
        J = 0 
        Icount = 0
        Num_Quarter_Atom = 0
        while I < self.NAtom:
            if (self.List_Num_Attached_Atom[I] == self.List_Num_Attached_Bond[I] and self.IAn[I] == Central_Atom_Type and self.List_Circle_Atom[I] == 0): # if the number of attached bonds equals to the number of attached atom, which means that all bonds are single bond
                for J in self.List_Surrounding_Atom[I][1:]:
#                    print(self.List_Surrounding_Atom[I][1:])
                    if ( self.IAn[J] not in self.Terminal_Atom): # this condition will be adjusted 
                        Icount = Icount + 1
                if (Icount == Num_Attaching_C_Atom and Icount == self.List_Num_Attached_Atom[I]):   # if there is no terminal atom and the attaching atom number equales the number of attaching unterminal_Atom 
                    Num_Quarter_Atom = Num_Quarter_Atom + 1
#                    print(I,Icount,Num_Attaching_C_Atom,"C_3")
                    self.Atom_Has_Been_Count.append(I)           #  adding the atom to the list of atom having been counted. 
                Icount = 0 
            I = I + 1
        return Num_Quarter_Atom
        
    


class Mutiple_Bond_Functional_Group(Group_Single_Bond):
    def __init__(self,NAtom,BO,IAn,ICoor):
        Group_Single_Bond.__init__(self,NAtom,BO,IAn,ICoor)
        self.List_Terminal_group_Central_Atom  =  Group_Single_Bond.Terminal_Group(self)
    

                                                                                                                                                                                                         #                       -
    def Double_Atom_Mutiple_Bond_Group(self,Central_Atom_Type,Attaching_Atom_Type,Bond_Type,Num_Bond_Side_atom_attaching,Num_Terminal_Atom_Attacing_The_Central_Atom = 0,Num_Terminal_Atom_Attacing_The_Side_Atom = 0):# eg (-)-c=N (6,7,3,3), this method is only suitable for double atom group like Cyano(-CN) 
        I = 0 
        J = 0
        Total_Num_Group = 0
        while I < self.NAtom:
#            print(self.IAn[I], self.List_Circle_Atom[I],"************")
            if (self.IAn[I] == Central_Atom_Type and self.List_Circle_Atom[I] == 0):
                while J < self.NAtom:
#                    if (self.BO[I][J] == Bond_Type and self.IAn[J] == Attaching_Atom_Type and self.List_Num_Attached_Bond[J] == Num_Bond_Side_atom_attaching and self.List_Num_of_terminal_atom_Surrounding[I] == Num_Terminal_Atom_Attacing_The_Central_Atom and self.List_Num_of_terminal_atom_Surrounding[J] == Num_Terminal_Atom_Attacing_The_Side_Atom): # if all qualifications match, the functional group exists like -CH=CH- (6,6,2,4,1,1)
                    if (self.BO[I][J] == Bond_Type and self.IAn[J] == Attaching_Atom_Type and self.List_Num_Attached_Bond[J] == Num_Bond_Side_atom_attaching and self.List_Circle_Atom[J] == 0): # if all qualifications match, the functional group exists like -CH=CH- (6,6,2,4,1,1)
#                        print(I,J,"CNCNCNCNNCN")
                        Total_Num_Group = Total_Num_Group + 1
                        self.Atom_Has_Been_Count.append(I)           #  adding the atom to the list of atom having been counted. 
                        self.Atom_Has_Been_Count.append(J)           #  adding the atom to the list of atom having been counted. 
                        # since this functional group has two atom 
                    J = J + 1
            I = I + 1 
            J = 0
        if (Central_Atom_Type == Attaching_Atom_Type and Num_Terminal_Atom_Attacing_The_Central_Atom == Num_Terminal_Atom_Attacing_The_Side_Atom):  # if Central_Atom and Attaching_Atom are the same and the surrounding condition is identical, if will count repeatly 
            Total_Num_Group = int(Total_Num_Group / 2)
        return Total_Num_Group
    



    def Detect_Nitro(self):  # This method is aim to detect the Nitro functional group -NO2 
        I = 0 
#        J1 = 0
#        J2 = 0 
        Num_Nitro = 0
        while I < self.NAtom:
            if (self.IAn[I] == self.N_Nuclei_Num and self.List_Circle_Atom[I] == 0 ):
                for J in self.List_Surrounding_Atom[I][1:] :  # for all atom who is attacing with N atom 
                    if (self.IAn[J] == self.O_Nuclei_Num and self.BO[I][J] == 2):
                        for K in self.List_Surrounding_Atom[I][1:] :  # for all atom who is attacing with N atom 
                            if (self.IAn[K] == self.O_Nuclei_Num and J != K and  self.BO[I][K] == 2):
#                                print(I,J,K)
                                Num_Nitro = Num_Nitro + 1    
                                self.Atom_Has_Been_Count.append(I)           #  adding the atom to the list of atom having been counted. 
                                self.Atom_Has_Been_Count.append(J)           #  adding the atom to the list of atom having been counted. 
                                self.Atom_Has_Been_Count.append(K)           #  adding the atom to the list of atom having been counted. 

            I = I + 1
        Num_Nitro = int(Num_Nitro)
        return Num_Nitro        
    

    
    def Detect_Thiocyanate(self) : # This method is designded for detecting the Thiocyanate (-)S-CN  
        I = 0 
        J = 0
        K = 0 
        L = 0
        CC_Single = 0
        Num_Thiocyanate = 0
        while I < self.NAtom:
            if (self.IAn[I] == self.N_Nuclei_Num and self.List_Circle_Atom[I] == 0): 
                while J < self.NAtom:
                    if (self.IAn[J] == self.C_Nuclei_Num and self.BO[I][J] == 3):
                        while K < self.NAtom:
                            if (self.IAn[K] ==  self.S_Nuclei_Num and self.BO[J][K] == 1 and self.List_Num_Attached_Bond[K] == 1 ):
                                Num_Thiocyanate =  Num_Thiocyanate  + 1
                                self.Atom_Has_Been_Count.append(I)           #  adding the atom to the list of atom having been counted.
                                self.Atom_Has_Been_Count.append(J)           #  adding the atom to the list of atom having been counted.
                                self.Atom_Has_Been_Count.append(K)           #  adding the atom to the list of atom having been counted.
                            K =  K + 1
                    J = J + 1
                    K = 0
            I = I + 1
            J = 0
            K = 0
        return Num_Thiocyanate
    


    def Detect_Nitrate(self):  # This method is designded for detecting the Nitrate(NO3-)
        I = 0
        J = 0 
        Icount = 0
        List_Tmp = []
        Num_Nitrate = 0
        while I < self.NAtom:
            if (self.IAn[I] == self.N_Nuclei_Num and self.List_Circle_Atom[I] == 0):
                List_Tmp.append(I)
                while J < self.NAtom:
                    if (self.IAn[J] == self.O_Nuclei_Num and self.BO[I][J] == 2):
                        Icount = Icount + 1
                        List_Tmp.append(J)
                    J = J + 1
                if (Icount == 3):
                    Num_Nitrate =  Num_Nitrate + 1
                    for M in List_Tmp:
                        self.Atom_Has_Been_Count.append(M)           #  adding the atom to the list of atom having been counted.
                J = 0
                Icount = 0 
                List_Tmp = []
            I = I + 1
        return Num_Nitrate 



    def Detect_Sulfate(self,Num_of_Double_Bond_O, Num_of_Single_Bond_O, Num_OH = 0):  # this method is aimed to detect the Sulfate (SO4 2-)
        I = 0
        J = 0
        K = 0
        Double_Bond_O = 0
        Single_Bond_O = 0
        Num_Sulfate   = 0
        List_Tmp = []
        OH = 0
#        CO_Single     = 0
        while I < self.NAtom:
            if (self.IAn[I] == self.S_Nuclei_Num and self.List_Circle_Atom[I] == 0):
#                print("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs")
                List_Tmp.append(I)
                while J < self.NAtom:
                    if (self.IAn[J] == self.O_Nuclei_Num and self.BO[I][J] == 2 ):
#                        print("S=O")
                        Double_Bond_O = Double_Bond_O + 1
                        List_Tmp.append(J)
#                    if (self.IAn[J] == self.O_Nuclei_Num and self.BO[I][J] == 1 and len(self.List_Surrounding_Atom[J]) == 2):
                    if (self.IAn[J] == self.O_Nuclei_Num and self.BO[I][J] == 1):
#                        while K < self.NAtom:
#                            if (self.IAn[K] == self.C_Nuclei_Num and self.BO[K][J] == 1 ):
#                                print("S-O-C",I,J,K,self.IAn[K])
                        for K in self.List_Surrounding_Atom[J][1:]:
                            if (self.IAn[K] == self.H_Nuclei_Num):
                                OH = OH + 1
                        Single_Bond_O = Single_Bond_O + 1
                        List_Tmp.append(J)
#                            K = K + 1
#                        K = 0
                    J = J + 1
                if (Double_Bond_O == Num_of_Double_Bond_O and Single_Bond_O == Num_of_Single_Bond_O and OH == Num_OH):
                    Num_Sulfate = Num_Sulfate + 1
                    for M in List_Tmp:
                        self.Atom_Has_Been_Count.append(M)           #  adding the atom to the list of atom having been counted.
                J = 0
                Double_Bond_O = 0
                Single_Bond_O = 0
            OH = 0
            List_Tmp = []
            I = I + 1
        return Num_Sulfate



    def Detect_phosphate(self,Num_of_Double_Bond_O, Num_of_Single_Bond_O, Num_OH = 0): # this function is designed to detect the phosphate group (PO4 3-)
        I = 0
        J = 0
        Double_Bond_O = 0
        Single_Bond_O = 0
        Num_phosphate   = 0
        OH        = 0
        List_Tmp = []
        while I < self.NAtom:
#            print("World")
            if (self.IAn[I] == self.P_Nuclei_Num and self.List_Circle_Atom[I] == 0):
                List_Tmp.append(I)
#                print("Hello")
                while J < self.NAtom:
                    if (self.IAn[J] == self.O_Nuclei_Num and self.BO[I][J] == 2 ):
                        Double_Bond_O = Double_Bond_O + 1
                        List_Tmp.append(J)
                    if (self.IAn[J] == self.O_Nuclei_Num and self.BO[I][J] == 1 ):
                        for K in self.List_Surrounding_Atom[J][1:]:
                            if (self.IAn[K] == self.H_Nuclei_Num):
                                OH = OH + 1
#                        print(I,J,OH)
                        Single_Bond_O = Single_Bond_O + 1
                        List_Tmp.append(J)
                    J = J + 1
                if (Double_Bond_O == Num_of_Double_Bond_O and Single_Bond_O == Num_of_Single_Bond_O and OH == Num_OH ):
                    Num_phosphate = Num_phosphate + 1
                    for M in List_Tmp:
                        self.Atom_Has_Been_Count.append(M)           #  adding the atom to the list of atom having been counted.
                J = 0
                Double_Bond_O = 0
                Single_Bond_O = 0
                OH        = 0 
                List_Tmp = []
            I = I + 1
        return Num_phosphate



    def Detect_Ester(self):
        I = 0
        J1 = 0   # atom connected with C
        J2 = 0   # atom connected with C(we want to find O atom) in double bond
        J3 = 0  # Atom in the Another side of oxyen
        Num_Ester = 0
        Num_Carboxyl = 0
        Num_Carboxylate_ion = 0
        CC_Single = 0
        CO_Double = 0
        CO_Single = 0
        OH_Single = 0  # another side of O is a C atom. 
        O_ = 0         #  Oxygen ion
        List_Tmp = []
        while I < self.NAtom :
            if (self.IAn[I] == 6 and self.List_Circle_Atom[I] == 0):
                List_Tmp.append(I)
                while J1 < self.NAtom:
                    if (self.BO[I][J1] == 2 and self.IAn[J1] == self.O_Nuclei_Num): #O atom connected with C in double bond
                        CO_Double = CO_Double + 1
                        List_Tmp.append(J1)
#                    if (self.BO[I][J1] == 1 and self.IAn[J1] != self.H_Nuclei_Num and self.IAn[J1] != self.O_Nuclei_Num): # C Atom is connected with C in single bond
#                        CC_Single = CC_Single + 1
                    J1 = J1 + 1
                
                while J2 < self.NAtom:
                    if (self.BO[I][J2] == 1 and self.IAn[J2] == self.O_Nuclei_Num): # O atom connected with C in single bondif (BO[I][J2] == 1 and IAn[J2] == 8): # O atom connected with C in single bond
                        CO_Single = CO_Single + 1
                        List_Tmp.append(J2)
                        if (self.List_Num_Attached_Atom[J2] == 1):
                            O_ = O_ + 1
#                            print("*********************")
                        else:
                            while J3 < self.NAtom:
                                if (self.BO[J2][J3] == 1 and J3 != I and self.IAn[J3] == self.H_Nuclei_Num): 
                                    OH_Single =  OH_Single + 1
                                if (self.BO[J2][J3] == 1 and J3 != I and self.IAn[J3] != self.H_Nuclei_Num): 
                                    OH_Single =  OH_Single + 0
                                J3 = J3 + 1
                            J3 = 0
                    J2 = J2 + 1
                J2 = 0
                J1 = 0
#            if (CO_Double == 1 and CO_Single == 1 and CC_Single == 1):  # if CO double bond and CO single bond exist simultaneously, Ester appears
            if (CO_Double == 1 and CO_Single == 1 ):  # if CO double bond and CO single bond exist simultaneously, Ester appears
                if (OH_Single == 1):
                    Num_Carboxyl = Num_Carboxyl + 1
#                    print("hello")
                    for M in List_Tmp:
                        self.Atom_Has_Been_Count.append(M)           #  adding the atom to the list of atom having been counted.
                elif (OH_Single == 0 and O_ == 0):
                    Num_Ester = Num_Ester + 1
#                    print ("Helooo")
                    for M in List_Tmp:
                        self.Atom_Has_Been_Count.append(M)           #  adding the atom to the list of atom having been counted.
                elif (O_ == 1 and OH_Single == 0):
                    Num_Carboxylate_ion = Num_Carboxylate_ion + 1
                    for M in List_Tmp:
                        self.Atom_Has_Been_Count.append(M)
            CO_Double = 0
            CO_Single = 0
            CC_Single = 0
            OH_Single = 0
            O_ = 0
            List_Tmp = []
            I = I + 1
        return Num_Ester,Num_Carboxyl,Num_Carboxylate_ion
                            

    def N_C_N(self):  #Cyanamide
        I = 0
        J = 0
#        K = 0
        Num_of_Cyanamide = 0
        Num_of_N3 = 0
        Num_of_N2 = 0
        Num_of_C2 = 0
        List_Tmp  = []
        while I < self.NAtom:
            if (self.IAn[I] == 6  and self.List_Num_Attached_Atom[I] == 2 and self.List_Num_Attached_Bond[I] and self.List_Circle_Atom[I] == 0):
#                print(self.List_Num_Attached_Bond[I])
#                print("CCCCCCCCCCC")
                Num_of_C2 = Num_of_C2 + 1
                List_Tmp.append(I)
                while J < self.NAtom:
                    if (self.IAn[J] == 7 and self.List_Circle_Atom[I] == 0 and self.BO[I][J] == 2 and self.List_Num_Attached_Bond[J] == 3 and self.List_Num_Attached_Atom[J] == 2):
                        Num_of_N3 = Num_of_N3 + 1
#                        print("N3")
                        List_Tmp.append(J)
                    if (self.IAn[J] == 7 and self.List_Circle_Atom[I] == 0 and self.BO[I][J] == 2 and self.List_Num_Attached_Bond[J] == 2 and self.List_Num_Attached_Atom[J] == 1):
                        Num_of_N2 = Num_of_N2 + 1
#                        print("N2")
                        List_Tmp.append(J)
                    J = J + 1
                J = 0
                if (Num_of_C2 == 1 and Num_of_N2 == 1 and Num_of_N3 == 1):
                    Num_of_Cyanamide = Num_of_Cyanamide + 1
                    for K in List_Tmp:
                        self.Atom_Has_Been_Count.append(K)
            Num_of_N3 = 0 
            Num_of_N2 = 0
            Num_of_C2 = 0
            List_Tmp  = []
            I = I + 1
        return Num_of_Cyanamide






if __name__ == "__main__":
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
    X=Group_Single_Bond(NAtom,BO,IAn,ICoor)
    print(X.IAn,X.List_Hybrid,X.List_Num_Attached_Atom,X.List_Num_Attached_Bond)
    _List_Terminal_group_Central_Atom = X.Terminal_Group()
    print(X.List_Circle_Atom)
    print("XXXXXXXXXXXXXXX")
    _CH3 = X.Detecting_Simple_Connecting_Group(6,1,3)
    _CH2 = X.Detecting_Simple_Connecting_Group(6,1,2)
    _PH3 = X.Detecting_Simple_Connecting_Group(15,1,3)
    _Number_Halogen_Atom = X.Halogen_Counting()
    __O_ = X.Quarter_Atom(8,2)
    print(_CH3,_CH2)
    print(X.Atom_Has_Been_Count)
#    print( X.List_Num_of_terminal_atom_Surrounding)
#    print(_CH3,_CH2,__O_)
    Y=Mutiple_Bond_Functional_Group(NAtom,BO,IAn,ICoor)
    _CO = Y.Double_Atom_Mutiple_Bond_Group(6,8,2,2)
    _CN = Y.Double_Atom_Mutiple_Bond_Group(6,7,3,3)
    _CHCH_ = Y.Double_Atom_Mutiple_Bond_Group(6,6,2,4,1,1)
    _NO2 = Y.Detect_Nitro()
    _NO3 = Y.Detect_Nitrate()
    _Sulfate = Y.Detect_Sulfate()
    _Num_Ester = Y.Detect_Ester()
    print(_CH3,_CO,_CN,_CHCH_,_NO2,_NO3,_Sulfate,_Num_Ester)
    _CH0 = X.Quarter_Atom(6,4)
    print(_CH0)
#    print(X.List_Hybrid,X.List_Surrounding_Atom,_List_Terminal_group_Central_Atom)




