import nupack as nu 
import numpy as np 
import pandas as pd 
import networkx as nx 
import re

from itertools import pairwise, tee

from .strand import Strand
from .model import Model

class Complex(object):
    
    def __init__(   self, 
                    model: Model, 
                    s1: Strand, 
                    s2: Strand, 
                    state,
                    structure=None,
                    dpxdist=None,
                    clean=False
                ):
        
        if type(model) != Model:
            raise TypeError("Model must be an instance of hdna.Model")
        self.model = model

        if type(s1) != Strand:
            raise TypeError("s1 must be an instance of hdna.Strand")
        if type(s2) != Strand:
            raise TypeError("s2 must be an instance of hdna.Strand")

        self.s1 = s1    #53
        self.s2 = s2    #35

        self.l1 = s1.length
        self.l2 = s2.length

        # Think of this as self.initialstructure
        # I didn't change it because of laziness 
        self.structure = structure

        self.dpxdist = dpxdist

        # self.getnupackproperties()
        if not clean:
            self.structureG(self.structure)        

        self.possible_states = [ None,
                            'singlestranded',
                            'duplex',
                            'zipping',
                            'on_nucleation',
                            'off_nucleation', 
                            'backfray',
                            'sliding' ]   

    
        if state not in self.possible_states:
            raise ValueError(f'state must be one among {self.possible_states} but you gave {state}')
        else: 
            self.state = state 
        
        if self.state == 'duplex':
            #TODO update this when considering mismatches 
            self.consecutive_nucleations = min(self.l1, self.l2)
            self.total_nucleations = self.consecutive_nucleations
            if self.structure == None:
                self.structure = '('*self.s1.length+'+'+')'*self.s2.length
            self.structureG()
            self.dpxdist = 0
        elif self.state == 'singlestranded':
            self.total_nucleations = 0
            self.consecutive_nucleations = 0
            self.structure = '.'*self.s1.length+'+'+'.'*self.s2.length
            self.G = 0
        else:
            self.total_nucleations = self.totbasepairs(self.structure)
            self.consecutive_nucleations = self.maxconsbp(self.structure)

        """ These variables are used to say if the current
        Coils object is duplexed, is offregister or is onregister
        to be comunicated to classes higher in the hierarchy"""

        # if self.state == 'on_nucleation': self.zipping_trajectory()
    
    def set_state(self, target):
        if target != self.state:
            self.state = target

        
    
######################################
##### Self-calculated Properties #####
######################################

    #TODO: ADD ABASIC DEFECTS (not supported by nupack) 
    # ab = {
    # 'GUG':  {'A':-4.3,'C':-3,'G':-5.7,'T':-4.5},
    # 'CUC':  {'A':-7.6,'C':-11.3,'G':-8.3,'T':-11.3},
    # }

    def totbasepairs(self, structure):
        nl = structure.count('(')
        nr = structure.count(')')
        if nl == nr: tot = nl
        else: raise BrokenPipeError(f'left and right nucleation should match: {structure}')
        return tot
    
    def maxconsbp(self, structure):
        l, r = structure.split('+')
        maxl = len(max(l.split('.')))
        maxr = len(max(r.split('.')))
        if maxl == maxr: cons = maxl
        else: raise BrokenPipeError(f'left and right nucleation should match: {structure}')
        return cons
        
    def inherit_zipping(self, listofzippings):
        self.zipping = listofzippings


#####################################
##### Nupack-dependant Methods  #####
#####################################

    def structureG(self, structure=None):
        """ Nupack related properties """
        if structure == None: structure = self.structure
        nuStrand1 = nu.Strand(self.s1.sequence, name = 'a')
        nuStrand2 = nu.Strand(self.s2.sequence, name = 'b') # An inversion here is needed because in this program strands are defined as 5-3 against 3-5 but in NUPACK all strands are defined 5-3 and the program takes care to turn them around and so on
        nuStructure = nu.Structure(structure)
        dG = nu.structure_energy(strands=[nuStrand1,nuStrand2], structure=nuStructure, model=self.model.nupack)
        self.G = dG
        return dG 


    def getnupackproperties(self):
            """ Nupack related properties """
            self.getnupackmodel()
            self.nuStrand1 = nu.Strand(self.s1.sequence, name = 'a')
            self.nuStrand2 = nu.Strand(self.s2.sequence, name = 'b') # An inversion here is needed because in this program strands are defined as 5-3 against 3-5 but in NUPACK all strands are defined 5-3 and the program takes care to turn them around and so on
            self.nuComplex = nu.Complex([self.nuStrand1,self.nuStrand2], name = 'c')
            j = nu.pfunc(self.nuComplex, self.model.nupack)
            self.duplexG = j[1]
            self.duplexZ = float(j[0])
            self.mfe = nu.mfe([self.nuStrand1, self.nuStrand2], model=self.model.nupack)


#############################################
##### General Structure Parsing Methods #####
#############################################
    
    """ This is a more general structure parsing routine 
        borrowed from the chamber class and needed here to
        parse zipping trajectories which for now don't account
        for wattson-crick base pairing.
        Since structures here have the normal dot-bracket notation
        this method doesn't work for ì structures. Maybe I will
        make a Structure class for handling all these structure
        related things """

    def parse_structure(self, structure, seq_a, seq_b):
        """
        Converts the internal general structures made with
        'ì' and '.' to dot-bracket structures considering the
        Wattson Crick pairings of bases
        """
        struct_a, struct_b = structure.split('+')
        seq_a = seq_a.sequence
        seq_b = seq_b.sequence 
        sx, dx = self.structurecut(seq_a, seq_b[::-1], struct_a, struct_b[::-1])
        patch_sx = ''.join(['\(' for i in range(len(sx))])
        patch_dx = ''.join(['\)' for i in range(len(dx))])
        out1 = re.sub(patch_sx,sx,struct_a)
        out2 = re.sub(patch_dx,dx,struct_b)
        structureout = out1+"+"+out2
        return structureout


    def structurecut(self, string1, string2, structure1, structure2):
        cut1 = ''.join([n1 for n1, s1 in zip(string1, structure1) if s1 == '('])
        cut2 = ''.join([n2 for n2, s2 in zip(string2, structure2) if s2 == ')'])
        sx = ''
        dx = ''
        # print(cut1, cut2)
        for n1, n2 in zip(cut1, cut2):
            if self._iswattsoncrick(n1, n2):
                sx += '('
                dx += ')'
            else: 
                sx += '.'
                dx += '.' 
        return sx, dx[::-1]

##########################
##### Helper Methods #####
##########################

    def _iswattsoncrick(self,n1,n2):
        """
        check for mismatches
        """
        if (n1 == 'A' and n2 == 'T') or (n1 == 'C' and n2 == 'G') or (n1 == 'T' and n2 == 'A') or (n1 == 'G' and n2 == 'C'):
            return True
        else:
            return False 

    def _get_mismatches(self):
        for nuc1, nuc2 in zip(self.s1.sequence, self.s2.sequence):
            self.mismatches.append(self._iswattsoncrick(nuc1, nuc2))

    def splitstructure(self):
        self.struct_a, self.struct_b = self.structure.split('+')
        return self.struct_a, self.struct_b

######################################
##### Methods Under Construction #####
######################################

    def correctedstrands(self):
        """ if given sequences have abasic defects this 
        method corrects the given Strands to non abasic 
        equivalent to pass them to nupack, the abasic
        penalty will be taken care of in the next method"""
        pass

    def abasicpenalty(self):
        """ Check the presence of an abasic defect
        and take the penalty out of the currently 
        calculated free energy """
        pass

######################
##### Properties #####
######################

    @property
    def sequences(self):
        return self.s1.sequence+"+"+self.s2.sequence



  
  