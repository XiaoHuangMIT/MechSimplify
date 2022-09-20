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



def find_central_coord_atom(molecule):
  
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
  
  
  
def find_coord_atom_distances(molecule,ccaid=None):
    
    #Analyzing distances between coordinating atoms
    #Analysis based on one ligand only since homoleptic
    #Inputs:
    #molecule: mol3D object (octahedral and homoleptic)
    #ccaid: one of the central coordinating atoms, default None for fac complexes
    #Outputs:
    #List of distances, and list of pairs: with the distance
    #If mer complex: two distances between central and two other coordinating atoms first
    
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

    return [dcc1,dcc2,dcc3],[(idx1,idx2),(idx1,idx3),(idx2,idx3)]

  
  
def find_cbond_length(molecule,ccaid=None):
  
    #Analyze bond length of coordination bonds
    #Analysis based on one ligand only since homoleptic
    #Inputs: 
    #molecule: mol3D object (octahedral and homoleptic)
    #ccaid: if mer complex
    #Outputs:
    #List of bond lengths and list of corresponding coordinating atoms
    #If mer: central coord atom first
    
    metal_id = molecule.findMetal()[0]
    if ccaid != None:
      lig1caid,lig2caid = ligand_breakdown(molecule)[2][0],ligand_breakdown(molecule)[2][1]
      if ccaid in lig1caid:
          idxs = lig1caid #indexes of coordinating atoms
      else:
          idxs = lig2caid
      idxs.remove(ccaid)
      idx1,idx2,idx3 = ccaid,idxs[0],idxs[1]
    else:
      idxs = ligand_breakdown(molcule)[2][0]
      idx1,idx2,idx3 = idxs[0],idxs[1],idxs[2]
      
    metal = molecule.getAtom(metal_id)
    atom1 = molecule.getAtom(idx1)
    atom2 = molecule.getAtom(idx2)
    atom3 = molecule.getAtom(idx3)
    dcm1 = metal.distance(atom1)
    dcm2 = metal.distance(atom2)
    dcm3 = metal.distance(atom3)
    
    return [dcm1,dcm2,dcm3],[idx1,idx2,idx3]

  
  
def find_cbond_angle(molecule,ccaid=None):

    #Analyze bond angles
    #Inputs:
    #molecule: mol3D object (octahedral, homoleptic)
    #ccaid: central coordinating atom (if mer symmstry)
    #Outputs:
    #List of bond angles, and list of pairs of corresponding coord atoms
    #If mer symmetry: two angles involving central coord atom first
    
    metal_id = molecule.findMetal()[0]
    if ccaid != None:
      lig1caid,lig2caid = ligand_breakdown(molecule)[2][0],ligand_breakdown(molecule)[2][1]
      if ccaid in lig1caid:
          idxs = lig1caid #indexes of coordinating atoms
      else:
          idxs = lig2caid
      idxs.remove(ccaid)
      idx1,idx2,idx3 = ccaid,idxs[0],idxs[1]
    else:
      idxs = ligand_breakdown(molecule)[2][0]
      idx1,idx2,idx3 = idxs[0],idxs[1],idxs[2]
      
     metal = molecule.getAtom(metal_id)
     atom1 = molecule.getAtom(idx1)
     atom2 = molecule.getAtom(idx2)
     atom3 = molecule.getAtom(idx3)
     v1 = np.array(metal.distancev(atom1))
     v2 = np.array(metal.distancev(atom2))
     v3 = np.array(metal.distancev(atom3))
     cmcangle1 = vecangle(v1,v2)
     cmcangle2 = vecangle(v1,v3)
     cmcangle3 = vecangle(v2,v3)
      
     return [cmcangle1,cmcangle2,cmcangle3],[(idx1,idx2),(idx1,idx3),(idx2,idx3)]
  
  
  
def 
  
  

