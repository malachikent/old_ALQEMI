"""
Unit and regression test for the alqemi package.
"""

# Import package, test suite, and other packages as needed
import sys
import pytest
from sympy import Dummy
from sympy.physics.secondquant import F, Fd, Commutator, KroneckerDelta
from sympy import sympify
from alqemi.core.tools import FermionAntiCommutator
from alqemi.algebra.rules import alqemi_reduce

def test_alqemi_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "alqemi" in sys.modules

# alqemi_reduce testing
def test_anticommutator_delta():
    """
    Tests Helgaker anticommutator identities
    Example: {p^+, q}
    """

    # 1. Define distinct dummy indices to prevent accidental trace collapses
    p, q = [Dummy(x) for x in "pq"]

    left_string = Fd(p)

    right_string = F(q)

    comm = FermionAntiCommutator(left_string,right_string)

    result = alqemi_reduce(comm)

    expected = KroneckerDelta(p,q)

    assert expected == result

def test_anticommutator_0():
    """
    Tests Helgaker anticommutator identities
    Example: {p, q}
    """

    # 1. Define distinct dummy indices to prevent accidental trace collapses
    p, q = [Dummy(x) for x in "pq"]

    left_string = F(p)

    right_string = F(q)

    comm = FermionAntiCommutator(left_string,right_string)

    result = alqemi_reduce(comm)

    expected = 0 

    assert expected == result
    
def test_1_2_commutator():
    """
    Tests the helgaker rules for [1, {1,1}] operator string.
    Example: [p^+, q^+ r]
    """

    # 1. Define distinct dummy indices to prevent accidental trace collapses
    p, q, r, s = [Dummy(x) for x in "pqrs"]

    # 2. Build the operator strings (length 2 = even parity)
    # Left: a^\dagger_p a_q
    left_string = Fd(p)

    # Right: a^\dagger_r a_s
    right_string = Fd(q) * F(r)

    # 3. Construct the unexpanded Commutator
    comm = Commutator(left_string, right_string)

    # 4. Execute your custom logic
    result = alqemi_reduce(comm)

    expected = -Fd(q) * KroneckerDelta(p,r) 

    # 6. Assert the result matches the expectation
    # See the note below regarding exact vs. algebraic equality
    assert result.doit() == expected

# @pytest.mark.even
def test_2_2_commutator():
    """
    Tests the helgaker_expansion rule for an [even, even] operator string.
    Example: [p^\dagger q, r^\dagger s] where both sides have 2 operators.
    """
    # 1. Define distinct dummy indices to prevent accidental trace collapses
    p, q, r, s = [Dummy(x) for x in "pqrs"]

    # 2. Build the operator strings (length 2 = even parity)
    # Left: a^\dagger_p a_q
    left_string = Fd(p) * F(q)
    
    # Right: a^\dagger_r a_s
    right_string = Fd(r) * F(s)

    # 3. Construct the unexpanded Commutator
    comm = Commutator(left_string, right_string)

    # 4. Execute your custom logic
    result = alqemi_reduce(comm)

    expected = KroneckerDelta(q,r)*Fd(p)*F(s) - KroneckerDelta(p,s)*Fd(r)*F(q) 

    # 6. Assert the result matches the expectation
    # See the note below regarding exact vs. algebraic equality
    assert result == expected

def test_2_4_commutator():
    """
    Tests the helgaker_expansion rule for an [even, even] operator string.
    Example: [p^+ q, r^+ s] where both sides have 2 operators.
    """
    # 1. Define distinct dummy indices to prevent accidental trace collapses
    p, q, r, s, m, n = [Dummy(x) for x in "pqrsmn"]

    # 2. Build the operator strings (length 2 = even parity)
    # Left: a^\dagger_p a_q
    left_string = Fd(p) * F(q)
    
    # Right: a^\dagger_r a_s
    right_string = Fd(r) * F(s) * Fd(m) * F(n)

    # 3. Construct the unexpanded Commutator
    comm = Commutator(left_string, right_string)

    # 4. Execute your custom logic
    result = alqemi_reduce(comm)
    print('RESULT:', result.expand())

    # 5. Construct the strictly expected SymPy tree for the even-even rule.
    expected = (
            KroneckerDelta(q,r)*Fd(p)*F(s)*Fd(m)*F(n) - 
            KroneckerDelta(p,s)*Fd(r)*F(q)*Fd(m)*F(n) + 
            KroneckerDelta(q,m)*Fd(r)*F(s)*Fd(p)*F(n) - 
            KroneckerDelta(p,n)*Fd(r)*F(s)*Fd(m)*F(q) 
            )
    # 6. Assert the result matches the expectation
    # See the note below regarding exact vs. algebraic equality
    assert expected.compare(result.doit()) == 0

def test_2_com_2_2_commutator():
    """
    Tests the helgaker expansion rule for an [even,[even,even]] operator string.
    Example: [p^+ q,[r^+ s, m^+ n]]
    """

    # 1. Define distinct dummy indices to prevent accidental trace collapses
    p, q, r, s, m, n = [Dummy(x) for x in "pqrsmn"]

    # 2. Build the operator strings (length 2 = even parity)
    left_string_1 = Fd(p) * F(q)

    left_string_2 = Fd(r) * F(s)

    right_string  = Fd(m) * F(n)

    comm = Commutator(left_string_1,Commutator(left_string_2, right_string))

    result = alqemi_reduce(comm).doit().expand()
    print('RESULT: ', result)

    expected = (
            KroneckerDelta(s,m)*KroneckerDelta(q,r)*Fd(p)*F(n) - 
            KroneckerDelta(s,m)*KroneckerDelta(p,n)*Fd(r)*F(q) - 
            KroneckerDelta(r,n)*KroneckerDelta(q,m)*Fd(p)*F(s) + 
            KroneckerDelta(r,n)*KroneckerDelta(p,s)*Fd(m)*F(q)
            )
    expected = expected.doit().expand()
    print('EXPECTED:', expected)

    assert expected.compare(result) == 0

def test_3_3_anticommutator():
    """
    Test the helgaker expansion rule for an [odd,odd]+ operator string.
    Example: [p^+ q r^+, s m^+ n]
    """
    # 1. Define distinct dummy indices to prevent accidental trace collapses
    p, q, r, s, m, n = [Dummy(x) for x in "pqrsmn"]

    left_string = Fd(p) * F(q) * Fd(r)

    right_string = F(s) * Fd(m) * F(n)

    comm = FermionAntiCommutator(left_string,right_string)

    result = alqemi_reduce(comm).doit().expand()
    print('RESULT: \n', result)
    expected = (
            KroneckerDelta(m, q) * F(s) * Fd(p) * Fd(r) * F(n) +
            KroneckerDelta(n, p) * F(s) * Fd(m) * F(q) * Fd(r) +
            KroneckerDelta(n, r) * F(s) * Fd(m) * Fd(p) * F(q) +
            KroneckerDelta(p, s) * F(q) * Fd(r) * Fd(m) * F(n) +
            KroneckerDelta(r, s) * Fd(p) * F(q) * Fd(m) * F(n)
        )
    expected = expected.doit().expand()
    print('EXPECTED: \n', expected)

    assert expected.compare(result) == 0


