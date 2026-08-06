from sympy import Sum, Add, Mul, Expr, sympify, S, Tuple, Function 
from sympy.core.sorting import default_sort_key
from sympy.physics.quantum import AntiCommutator, Commutator
from sympy.physics.secondquant import CreateFermion, AnnihilateFermion, F, Fd, evaluate_deltas
from sympy import KroneckerDelta
from sympy.combinatorics import Permutation

@property
def is_left_killer(self):

    if isinstance(self, CreateFermion) and self.is_above_fermi:
        return True
    
    if isinstance(self,AnnihilateFermion) and self.is_below_fermi:
        return True

    else:
        return False

@property 
def is_right_killer(self):

    if isinstance(self, CreateFermion) and self.is_below_fermi:
        return True
    
    if isinstance(self,AnnihilateFermion) and self.is_above_fermi:
        return True

    else:
        return False
AnnihilateFermion.is_left_killer = is_left_killer
AnnihilateFermion.is_right_killer = is_right_killer

CreateFermion.is_left_killer = is_left_killer
CreateFermion.is_right_killer = is_right_killer

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

class DoubleBar(Function):
    """
    Symbolic representation of an antisymmetrized two-electron integral: <pq || rs>
    """
    
    @classmethod
    def eval(cls, p, q, r, s):
        # 1. Pauli Exclusion: If any two indices in the same bra/ket match, it is strictly 0
        if p == q or r == s:
            return S.Zero
            
        sign = 1
        p1, q1, r1, s1 = p, q, r, s
        
        # 2. Antisymmetry in the bra: <qp || rs> = -<pq || rs>
        # We sort them lexicographically so 'p' always comes before 'q'
        if default_sort_key(p1) > default_sort_key(q1):
            p1, q1 = q1, p1
            sign *= -1
            
        # 3. Antisymmetry in the ket: <pq || sr> = -<pq || rs>
        if default_sort_key(r1) > default_sort_key(s1):
            r1, s1 = s1, r1
            sign *= -1
            
        # 4. Hermitian Symmetry (for real orbitals): <rs || pq> = <pq || rs>
        # We compare the bra pair to the ket pair and sort them
        pair1 = (default_sort_key(p1), default_sort_key(q1))
        pair2 = (default_sort_key(r1), default_sort_key(s1))
        
        if pair1 > pair2:
            p1, q1, r1, s1 = r1, s1, p1, q1
            # Swapping the bra and ket does not induce a sign change
            
        # 5. The Output Decision
        # If the indices had to be rearranged, return the newly sorted object multiplied by the sign
        if (p1, q1, r1, s1) != (p, q, r, s):
            return sign * cls(p1, q1, r1, s1)
            
        # If the indices were already perfectly sorted, return None. 
        # This tells SymPy: "Everything is fine, instantiate the function normally."
        return None

    # --- Pretty Printing Methods ---
    
    def _sympystr(self, printer):
        """Terminal string representation: <pq||rs>"""
        p, q, r, s = [printer._print(arg) for arg in self.args]
        return f"<{p}{q}||{r}{s}>"
        
    def _latex(self, printer):
        """Jupyter Notebook / LaTeX representation: \langle pq || rs \rangle"""
        p, q, r, s = [printer._print(arg) for arg in self.args]
        return rf"\langle {p} {q} || {r} {s} \rangle"
