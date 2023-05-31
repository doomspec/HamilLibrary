from hamil_lib_helper import *
import tequila as tq

if verify_hash(__file__):
    exit(0)

basis_set = 'sto-6g'
molecule = tq.quantumchemistry.Molecule(geometry="./geometry/H2O.xyz", basis_set=basis_set)

save_fermion_hamil_for_mol(molecule, basis_set, __file__)
