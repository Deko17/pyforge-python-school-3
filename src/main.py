from rdkit import Chem

class Molecule:
    def __init__(self):
        self.molecules = {}
    
    def add_molecules(self,identifier,smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            self.molecules[identifier] = smiles
        else:
            raise ValueError("Invalid string")
        
    def get_molecules(self, identifier):
        return self.molecules.get(identifier, None)
    
    def update_molecules(self, identifier, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            self.molecules[identifier] = smiles
        else:
            raise ValueError ("Invaid string")
        
    def delete_molecules(self, identifier):
        if identifier in self.molecules:
            del self.molecules[identifier]
        else:
            raise ValueError ("Identifier invalid")
    
    def list_molecules(self):
        return self.molecules
    
    def substructure_search(smiles_list, sub_smiles):
        sub_mol = Chem.MolFromSmiles(sub_smiles) #for conversion substructure SMILES to RDKIT molecujle
        matches = []        #empty list for storage of molecules
        for smiles in smiles_list:  #converts each smilel to rdkit mol and if it is created it checks if it is valid structure
            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.HasSubstructMatch(sub_mol):      
            #adds molecules that are matched to list
                matches.append(smiles)
        return matches
    

mol = Molecule()

mol.add_molecules("mol1", "CCO")
mol.add_molecules("mol2", "c1ccccc1")
mol.add_molecules("mol3", "CC(=O)O")

print("Get function:", mol.get_molecules("mol1"))

print("List function:", mol.list_molecules())

mol.delete_molecules("mol2")

mol.update_molecules("mol2", "CCCOCCC")

substructure = "c1ccccc1"
result = mol.substructure_search(substructure)
print(result) 
