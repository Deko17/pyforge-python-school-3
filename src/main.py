import logging
import json
from fastapi import FastAPI, HTTPException, Query
from rdkit import Chem
import redis

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI()

redis_client = redis.Redis(host='redis', port=6379, db=0)

CACHE_EXPIRATION = 60  

class MoleculeManager:
    def __init__(self):
        self.molecules = {}

    def add_molecule(self, identifier, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            self.molecules[identifier] = smiles
            logger.info(f"Added molecule with identifier: {identifier}")
            redis_client.delete('molecules_list')
        else:
            logger.error(f"Failed to add molecule with identifier: {identifier}. Invalid SMILES string.")
            raise ValueError("Invalid SMILES string")

    def get_molecule(self, identifier):
        logger.info(f"Fetching molecule with identifier: {identifier}")
        return self.molecules.get(identifier, None)

    def update_molecule(self, identifier, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            self.molecules[identifier] = smiles
            logger.info(f"Updated molecule with identifier: {identifier}")
            redis_client.delete('molecules_list')
        else:
            logger.error(f"Failed to update molecule with identifier: {identifier}. Invalid SMILES string.")
            raise ValueError("Invalid SMILES string")

    def delete_molecule(self, identifier):
        if identifier in self.molecules:
            del self.molecules[identifier]
            logger.info(f"Deleted molecule with identifier: {identifier}")
            redis_client.delete('molecules_list')
        else:
            logger.error(f"Failed to delete molecule with identifier: {identifier}. Not found.")
            raise KeyError("Identifier not found")

    def list_molecules(self, limit=None):
        logger.info(f"Listing molecules with limit: {limit}")
        
        cached_list = redis_client.get('molecules_list')
        if cached_list:
            logger.info("Returning cached molecules list")
            molecules = json.loads(cached_list)
        else:
            molecules = list(self.molecules.items())
            redis_client.setex('molecules_list', CACHE_EXPIRATION, json.dumps(molecules))
            logger.info("Cached molecules list")

        def molecule_iterator():
            for index, molecule in enumerate(molecules):
                if limit is not None and index >= limit:
                    break
                yield molecule
        return molecule_iterator()

    def substructure_search(self, sub_smiles):
        logger.info(f"Performing substructure search for: {sub_smiles}")

        cache_key = f"substructure_search:{sub_smiles}"
        cached_result = redis_client.get(cache_key)
        if cached_result:
            logger.info("Returning cached substructure search result")
            return json.loads(cached_result)

        sub_mol = Chem.MolFromSmiles(sub_smiles)
        if not sub_mol:
            logger.error(f"Substructure search failed. Invalid substructure SMILES string: {sub_smiles}")
            raise ValueError("Invalid substructure SMILES string")

        matches = []
        for identifier, smiles in self.molecules.items():
            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.HasSubstructMatch(sub_mol):
                matches.append((identifier, smiles))

        redis_client.setex(cache_key, CACHE_EXPIRATION, json.dumps(matches))
        logger.info(f"Substructure search completed with {len(matches)} matches. Cached result.")
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
def list_molecules(limit: int = Query(None, description="Limit the number of molecules returned")):
    return list(manager.list_molecules(limit))

@app.post("/molecules/search/")
def substructure_search(sub_smiles: str):
    try:
        results = manager.substructure_search(sub_smiles)
        return results
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
