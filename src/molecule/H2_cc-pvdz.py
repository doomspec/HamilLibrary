from hamil_lib_helper import *
import tequila as tq

basis_set = 'cc-pvdz'
molecule = tq.quantumchemistry.Molecule(geometry="./geometry/H2.xyz", basis_set=basis_set)

save_fermion_hamil(molecule, basis_set, __file__)

