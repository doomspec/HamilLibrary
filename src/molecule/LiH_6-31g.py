from hamil_lib_helper import *
import tequila as tq

basis_set = '6-31g'
molecule = tq.quantumchemistry.Molecule(geometry="./geometry/LiH.xyz", basis_set=basis_set)

save_fermion_hamil(molecule, basis_set, __file__)

