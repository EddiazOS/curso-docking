from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def analizar_molecula(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        # Calcular propiedades
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        h_donors = Lipinski.NumHDonors(mol)
        h_acceptors = Lipinski.NumHAcceptors(mol)
        rotatableble_bonds = Lipinski.NumRotatableBonds(mol)

        print("Peso Molecular: ", mw)
        print("Log_P: ", logp)
        print("Enlace Rotable: ", rotatableble_bonds)
        print("H donables: ", h_donors)
        print("H aceptables: ", h_acceptors,sep="\t")
    else:
        print("SMILE INVALIDO")

# Ejemplo: Aspirina
aspirina = "CC(=O)OC1=CC=CC=C1C(=O)O"
analizar_molecula(smiles=aspirina)