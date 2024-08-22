from fastapi import FastAPI, HTTPException
from rdkit import Chem

app = FastAPI()

class MoleculeManager:
    def __init__(self):
        self.molecules = {}

    def add_molecule(self, identifier, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            self.molecules[identifier] = smiles
        else:
            raise ValueError("Invalid SMILES string")

    def get_molecule(self, identifier):
        return self.molecules.get(identifier, None)

    def update_molecule(self, identifier, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            self.molecules[identifier] = smiles
        else:
            raise ValueError("Invalid SMILES string")

    def delete_molecule(self, identifier):
        if identifier in self.molecules:
            del self.molecules[identifier]
        else:
            raise KeyError("Identifier not found")

    def list_molecules(self):
        return self.molecules

    def substructure_search(self, sub_smiles):
        sub_mol = Chem.MolFromSmiles(sub_smiles)
        if not sub_mol:
            raise ValueError("Invalid substructure SMILES string")
        matches = []
        for identifier, smiles in self.molecules.items():
            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.HasSubstructMatch(sub_mol):
                matches.append((identifier, smiles))
        return matches

manager = MoleculeManager()

@app.post("/molecules/")
def add_molecule(identifier: str, smiles: str):
    try:
        manager.add_molecule(identifier, smiles)
        return {"message": "Molecule added successfully"}
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.get("/molecules/{identifier}")
def get_molecule(identifier: str):
    molecule = manager.get_molecule(identifier)
    if molecule:
        return {"identifier": identifier, "smiles": molecule}
    else:
        raise HTTPException(status_code=404, detail="Molecule not found")

@app.put("/molecules/{identifier}")
def update_molecule(identifier: str, smiles: str):
    try:
        manager.update_molecule(identifier, smiles)
        return {"message": "Molecule updated successfully"}
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.delete("/molecules/{identifier}")
def delete_molecule(identifier: str):
    try:
        manager.delete_molecule(identifier)
        return {"message": "Molecule deleted successfully"}
    except KeyError as e:
        raise HTTPException(status_code=404, detail=str(e))

@app.get("/molecules/")
def list_molecules():
    return manager.list_molecules()

@app.post("/molecules/search/")
def substructure_search(sub_smiles: str):
    try:
        results = manager.substructure_search(sub_smiles)
        return results
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
