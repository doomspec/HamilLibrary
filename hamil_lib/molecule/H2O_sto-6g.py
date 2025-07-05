from hamil_lib.hamil_lib_helper import *
import tequila as tq

basis_set = 'sto-6g'
molecule = tq.quantumchemistry.Molecule(geometry="./geometry/H2O.xyz", basis_set=basis_set)

save_fermion_hamil(molecule, basis_set, __file__)
