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
  
