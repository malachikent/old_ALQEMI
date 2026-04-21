"""
Unit and regression test for the alqemi package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import alqemi


def test_alqemi_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "alqemi" in sys.modules
