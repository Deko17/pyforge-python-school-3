from app.celery_maker import celery
from rdkit import Chem

@celery.task
def substructure_search_task(sub_smiles, molecules):
    sub_mol = Chem.MolFromSmiles(sub_smiles)
    if not sub_mol:
        raise ValueError("Invalid substructure SMILES string")

    matches = []
    for identifier, smiles in molecules.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.HasSubstructMatch(sub_mol):
            matches.append((identifier, smiles))

    return matches
