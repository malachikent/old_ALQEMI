from sympy import Sum, Add, Mul, Expr, sympify, S, Tuple 
from sympy.core.sorting import default_sort_key
from sympy.physics.quantum import AntiCommutator, Commutator
from sympy.physics.secondquant import F, Fd, wicks, evaluate_deltas
from sympy import KroneckerDelta
from sympy.combinatorics import Permutation

class StrictCommutator(Commutator):
    """ Bypasses SymPy's native Commutator auto-evaluating to 0. """
    is_commutative = False
    
    def __new__(cls, a, b):
        # By directly returning Expr.__new__, we skip the 0-check completely!
        return Expr.__new__(cls, sympify(a), sympify(b))

    def _eval_expand_commutator(self, **hints):
        # Teaches SymPy how to physically tear down the bracket!
        A, B = self.args
        return A * B - B * A

class FermionAntiCommutator(AntiCommutator):
    """
    Custom AntiCommutator that evaluates fundamental fermion identities 
    before falling back to standard algebraic expansion.
    """
    def __new__(cls, a, b):
        a = sympify(a)
        b = sympify(b)
        # Bypass SymPy's sort_key canonicalization!
        return Expr.__new__(cls, a, b)
    
    def doit(self, **hints):
        A, B = self.args

        #[A+,A+]+ = 0 or [A,A] = 0
        if A.func == B.func and A.func in (F, Fd):
            return 0
        
        #Anticommutator Identities of all permutations of elementary operators and orbitals(general, occupied, virtual)
        if A.func != B.func and A.func in (F, Fd):
            p = A.args[0]
            q = B.args[0]

            p_virt = p.assumptions0.get('above_fermi', False)
            p_occ = p.assumptions0.get('below_fermi', False)
            
            q_virt = q.assumptions0.get('above_fermi', False)
            q_occ = q.assumptions0.get('below_fermi', False)
            if (p_virt and q_occ) or \
               (p_occ and q_virt):
               return 0
            
            return KroneckerDelta(p,q, evaluate=False)
        # Fallback: If it's a string of operators, do the standard AB + BA
        return super().doit(**hints)

class DoubleBarIntegral(Expr):
    """
    Represents the antisymmetrized two-electron integral <pq||rs>.
    Self-canonicalizes upon creation and enforces Pauli exclusion.
    """
    is_commutative = True 

    def __new__(cls, *args):
        # Safely handle SymPy's internal argument flattening!
        if len(args) == 2:
            bra_indices, ket_indices = args
        elif len(args) == 4:
            bra_indices = args[:2]
            ket_indices = args[2:]
        else:
            raise ValueError(f"DoubleBarIntegral expects 2 Tuples or 4 Indices. Got: {args}")

        bra = tuple(sympify(i) for i in bra_indices)
        ket = tuple(sympify(i) for i in ket_indices)

        # 1. Pauli Exclusion Check
        if len(set(bra)) < len(bra) or len(set(ket)) < len(ket):
            return S.Zero

        # 2. Sort using SymPy's native sorting function
        bra_sorted = tuple(sorted(bra, key=default_sort_key))
        ket_sorted = tuple(sorted(ket, key=default_sort_key))

        # 3. Determine the parity of the permutations
        bra_perm = Permutation([bra.index(x) for x in bra_sorted])
        ket_perm = Permutation([ket.index(x) for x in ket_sorted])
        
        sign = S.NegativeOne**(bra_perm.parity() + ket_perm.parity())

        # 4. Create the strictly ordered object using SymPy Tuples
        obj = Expr.__new__(cls, Tuple(*bra_sorted), Tuple(*ket_sorted))
        
        return sign * obj

    def _latex(self, printer):
        bra_str = "".join(printer._print(i) for i in self.args[0])
        ket_str = "".join(printer._print(i) for i in self.args[1])
        return rf"\langle {bra_str} || {ket_str} \rangle"

