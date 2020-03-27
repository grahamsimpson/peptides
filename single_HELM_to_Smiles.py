# SMILES_From_HELM.py
from rdkit import Chem
from rdkit.Chem import AllChem

HELM = "PEPTIDE1{[ac].I.C.E.C.R.E.A.M.A.A.K.[am]}|PEPTIDE2{[ac].I.C.E.C.R.E.A.M.D.D}$PEPTIDE2,PEPTIDE1,11:R2-12:R3$$$"

mol = Chem.MolFromHELM(HELM)
print(mol)

if mol != None:
    smiles = Chem.MolToSmiles(mol)
    print(smiles)
else:
    print("error, MolFromHELM gave "None" - please check your HELM")
