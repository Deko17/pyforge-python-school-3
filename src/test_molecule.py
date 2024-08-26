import pytest
import logging
from main import main  # Replace `main` with your actual module name if it's different

# Configure logging for tests
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_add_molecule():
    logger.info("Starting test_add_molecule")
    manager = main()
    manager.add_molecule("mol1", "CCO")
    assert manager.get_molecule("mol1") == "CCO"
    logger.info("test_add_molecule passed")

def test_invalid_smiles():
    logger.info("Starting test_invalid_smiles")
    manager = main()
    with pytest.raises(ValueError):
        manager.add_molecule("mol2", "invalid_smiles")
    logger.info("test_invalid_smiles passed")

def test_empty_database():
    logger.info("Starting test_empty_database")
    manager = main()
    assert manager.list_molecules() == []
    logger.info("test_empty_database passed")

def test_get_nonexistent_molecule():
    logger.info("Starting test_get_nonexistent_molecule")
    manager = main()
    with pytest.raises(KeyError):
        manager.get_molecule("nonexistent_mol")
    logger.info("test_get_nonexistent_molecule passed")

def test_remove_molecule():
    logger.info("Starting test_remove_molecule")
    manager = main()
    manager.add_molecule("mol1", "CCO")
    manager.remove_molecule("mol1")
    with pytest.raises(KeyError):
        manager.get_molecule("mol1")
    logger.info("test_remove_molecule passed")

# def test_duplicate_molecule():
#     logger.info("Starting test_duplicate_molecule")
#     manager = main()
#     manager.add_molecule("mol1", "CCO")
#     with pytest.raises(ValueError):
#         manager.add_molecule("mol1", "CCO")
#     logger.info("test_duplicate_molecule passed")

def test_update_molecule():
    logger.info("Starting test_update_molecule")
    manager = main()
    manager.add_molecule("mol1", "CCO")
    manager.update_molecule("mol1", "CCC")
    assert manager.get_molecule("mol1") == "CCC"
    logger.info("test_update_molecule passed")
