import pytest
from main import main  

def test_add_molecule():
    manager = main()
    manager.add_molecule("mol1", "CCO")
    assert manager.get_molecule("mol1") == "CCO"

def test_invalid_smiles():
    manager = main()
    with pytest.raises(ValueError):
        manager.add_molecule("mol2", "invalid_smiles")
