import sys
from pathlib import Path

# Makes `from analysis import ...` and `from load_data import ...` work from
# inside tests/, without needing the project installed as a package.
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import sqlite3

import pytest


@pytest.fixture(scope="session")
def ensure_db() -> Path:
    """Builds cell_counts.db from the real cell-count.csv checked into the
    repo, once per test session, using the actual load_data.py pipeline
    (not a separate test-only copy of the loading logic)."""
    from load_data import build_database
    build_database()
    return ROOT / "cell_counts.db"


@pytest.fixture(scope="session")
def conn(ensure_db):
    """A real connection to the built database, shared across tests in a
    session for speed. Tests should not mutate the database."""
    connection = sqlite3.connect(ensure_db)
    connection.execute("PRAGMA foreign_keys = ON")
    yield connection
    connection.close()
