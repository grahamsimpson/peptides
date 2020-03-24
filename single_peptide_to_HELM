from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolToSequence
from numpy import unique, argmin
from numpy import abs as np_abs
from Rules import pkatable, condensedtable

input = "(Ac-ICECREAMAA)(Ac-ICECREAMDD)K-NH2"

# puts dots in between all characters if number of hyphens =2; if not replaces brackets by -
a = ".".join(input) if input.count(
    "-") == 2 else input.replace("(", "-").replace(")", "-")
# for linear peptides replaces --Ac- by -
b = a.replace('--Ac-', '') if a.count("-") == 7 else a
# for branched peptides inverts the second sequence for first HELM attempt
c = ".".join(b[0: 14] + b[14:26][::-1] + b[26:]
             ) if b.count("-") >= 4 else b
# for branched peptides inserts .X. after first sequence for second HELM attempt
d = ".".join(b[0:14] + ".X." + b[14:]) if b.count("-") >= 4 else b

print("this is c ", c, "\n")
print("this is d ", d, "\n")

HELM_string1 = c.replace("-.A.c.", "PEPTIDE1{[am]}|PEPTIDE2{[ac]}|PEPTIDE3{[ac]").replace(".r", ".[dR]").replace(".-", "").replace("-", "").replace(".f", ".[dF]").replace(".v", ".[dV]").replace(".[.N.l.e.]", ".[Nle]").replace(".a", ".[dA]").replace(".w", ".[dW]").replace(".i", ".[dI]").replace(".p", ".[dP]").replace(".s", ".[dS]").replace("['A.c.-.", "PEPTIDE1{[ac].").replace("A.c.-.", "PEPTIDE1{[ac].").replace(".N.H.2", "}$PEPTIDE2,PEPTIDE3,1:R2-22:R2|PEPTIDE1,PEPTIDE3,1:R1-12:R3$$$V2.0").replace(
    ".[.H.A.[dR].g.]", ".[Har]").replace(".[.b.A.l.a.]", ".[Bal]") if c.count("-") >= 4 else c.replace(".r", ".[dR]").replace(".f", ".[dF]").replace(".v", ".[dV]").replace(".[.N.l.e.]", ".[Nle]").replace(".a", ".[dA]").replace(".w", ".[dW]").replace(".i", ".[dI]").replace(".p", ".[dP]").replace(".s", ".[dS]").replace("['A.c.-.", "PEPTIDE1{[ac].").replace("A.c.-.", "PEPTIDE1{[ac].").replace("-.N.H.2", "[am]}$$$$V2.0").replace(".[.H.A.[dR].g.]", ".[Har]").replace(".[.b.A.l.[dA].]", ".[Bal]")

HELM_string2 = d.replace("-.A.c.", "PEPTIDE1{[ac].").replace(".r", ".[dR]").replace(".-", "").replace("-", "").replace(".f", ".[dF]").replace(".v", ".[dV]").replace(".[.N.l.e.]", ".[Nle]").replace(".a", ".[dA]").replace(".w", ".[dW]").replace(".i", ".[dI]").replace(".p", ".[dP]").replace(".s", ".[dS]").replace("['A.c.-.", "PEPTIDE1{[ac].").replace("A.c.-.", "PEPTIDE1{[ac].").replace("...X...", "}|PEPTIDE2{[ac].").replace(".N.H.2", ".[am]}$PEPTIDE2,PEPTIDE1,12:R3-11:R2$$$V2.0").replace(
    ".[.H.A.[dR].g.]", ".[Har]").replace(".[.b.A.l.a.]", ".[Bal]") if d.count("-") >= 4 else d.replace(".r", ".[dR]").replace(".f", ".[dF]").replace(".v", ".[dV]").replace(".[.N.l.e.]", ".[Nle]").replace(".a", ".[dA]").replace(".w", ".[dW]").replace(".i", ".[dI]").replace(".p", ".[dP]").replace(".s", ".[dS]").replace("['A.c.-.", "PEPTIDE1{[ac].").replace("A.c.-.", "PEPTIDE1{[ac].").replace("-.N.H.2", "[am]}$$$$V2.0").replace(".[.H.A.[dR].g.]", ".[Har]").replace(".[.b.A.l.[dA].]", ".[Bal]")

print("This is your HELM1 ", HELM_string1, "\n")
print("This is your HELM2 ", HELM_string2, "\n")

mol = Chem.MolFromHELM(HELM_string1)
Sequence = Chem.MolToSequence(mol)
SMILES = Chem.MolToSmiles(mol)
HELM = Chem.MolToHELM(mol)

print("This is your input ", input, "\n")
print("This is your mol ", mol, "\n")
print("This is your sequence ", sequence, "\n")
print("This is your smiles ", SMILES, "\n")
print("This is your HELM ", HELM_string, "\n")
