import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from molSimplifyAD.utils.pymongo_tools import connect2db, insert, count_find
from molSimplify.Classes.atom3D import atom3D
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.ligand import ligand_breakdown 
from molSimplify.Scripts.geometry import vecangle
import networkx as nx



def find_central_atom(molecule):
  
    #Returns the central coordinating atom of two ligands
    #Inputs:
    #molecule: mol3D object
    
    try:
      metal_id = molecule.findMetal()[0]
      lig1caid,lig2caid = ligand_breakdown(molecule)[2][0],ligand_breakdown(molecule)[2][1]
      lig1ca3d,lig2ca3d = [],[]
      for i in np.arange(3):
          lig1ca3d.append(molecule.getAtom(lig1caid[i]))
          lig2ca3d.append(molecule.getAtom(lig2caid[i]))
      lig1cadists,lig2cadists = [],[]
      lig1cadists.append(lig1ca3d[0].distance(lig1ca3d[1])+lig1ca3d[0].distance(lig1ca3d[2]))
      lig1cadists.append(lig1ca3d[1].distance(lig1ca3d[0])+lig1ca3d[1].distance(lig1ca3d[2]))
      lig1cadists.append(lig1ca3d[2].distance(lig1ca3d[0])+lig1ca3d[2].distance(lig1ca3d[1]))
      lig2cadists.append(lig2ca3d[0].distance(lig2ca3d[1])+lig2ca3d[0].distance(lig2ca3d[2]))
      lig2cadists.append(lig2ca3d[1].distance(lig2ca3d[0])+lig2ca3d[1].distance(lig2ca3d[2]))
      lig2cadists.append(lig2ca3d[2].distance(lig2ca3d[0])+lig2ca3d[2].distance(lig2ca3d[1]))
      ccaid1 = lig1caid[lig1cadists.index(min(lig1cadists))]#idx of 'central' atom in first and second ligand
      ccaid2 = lig2caid[lig2cadists.index(min(lig2cadists))]
    except:
      ccaid1, ccaid2 = 'Not found', 'Not found'
      
    return ccaid1,ccaid2
  
  
  
def find_dcc(molecule,ccaid=None):
    
    #Analyzing distances between coordinating atoms
    #Returns a list of list, element in format [idx1,idx2,distance]
    #Analysis based on one ligand only since homoleptic
    #If mer complex: first two elements correspond to distance between central and two other coordinating atoms
    #Inputs:
    #molecule: mol3D object
    #ccaid one of the central coordinating atoms, default None for fac complexes
    
    if ccaid != None: #Mer complexes
      metal_id = molecule.findMetal()[0]
      lig1caid,lig2caid = ligand_breakdown(molecule)[2][0],ligand_breakdown(molecule)[2][1]
      if ccaid in lig1caid:
        idxs = lig1caid #indexes of coordinating atoms
      else:
        idxs = lig2caid
      idxs.remove(ccaid)
      idx1, idx2, idx3 = ccaid,idxs[0],idxs[1] #central coordinating atom first
    else: #Fac complexes
      idxs = ligand_breakdown(molcule)[2][0]
      idx1, idx2, idx3 = idxs[0], idxs[1], idxs[2]
      
    atom1 = molecule.getAtom(idx1)
    atom2 = molecule.getAtom(idx2)
    atom3 = molecule.getAtom(idx3)
    dcc1 = atom1.distance(atom2)
    dcc2 = atom1.distance(atom3)
    dcc3 = atom2.distance(atom3)

    return [[idx1,idx2,dcc1],[idx1,idx3,dcc2],[idx2,idx3,dcc3]]

def 
