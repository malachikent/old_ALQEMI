from sympy import Symbol, Dummy, oo, Sum, Mul, Add, Expr, sympify, S, Rational, factorial, Tuple
from sympy.core.sorting import default_sort_key
from sympy.physics.secondquant import F, Fd
from sympy.combinatorics import Permutation
from core.tools import DoubleBarIntegral

class ClusterAmplitude(Expr):
    """
    Represents the antisymmetric cluster amplitude tensor t_{ij...}^{ab...}.
    Automatically canonicalizes indices and enforces Pauli exclusion.
    """
    is_commutative = True 
    
    def __new__(cls, label, holes, parts):
        holes = tuple(sympify(h) for h in holes)
        parts = tuple(sympify(p) for p in parts)
        
        # 1. Pauli Exclusion Check: If an index is repeated, the amplitude is exactly zero
        if len(set(holes)) < len(holes) or len(set(parts)) < len(parts):
            return S.Zero
            
        # 2. Sort using native SymPy logic for canonical memory mapping
        h_sorted = tuple(sorted(holes, key=default_sort_key))
        p_sorted = tuple(sorted(parts, key=default_sort_key))
        
        # 3. Determine the permutation parity
        h_perm = Permutation([holes.index(h) for h in h_sorted])
        p_perm = Permutation([parts.index(p) for p in p_sorted])
        
        sign = S.NegativeOne**(h_perm.parity() + p_perm.parity())
        
        # 4. Construct and return with the proper fermion phase
        obj = Expr.__new__(cls, Symbol(str(label)), Tuple(*h_sorted), Tuple(*p_sorted))
        return sign * obj

    def _latex(self, printer):
        # Renders beautifully in Jupyter: t_{ij}^{ab}
        h_str = "".join(printer._print(i) for i in self.args[1])
        p_str = "".join(printer._print(a) for a in self.args[2])
        return f"{self.args[0]}_{{{h_str}}}^{{{p_str}}}"


class ClusterOperator(Expr):
    """
    Abstract representation of a full cluster operator T_n.
    Generates unique concrete dummy indices only when evaluated.
    """
    is_commutative = False
    
    def __new__(cls, label, rank):
        obj = Expr.__new__(cls, Symbol(str(label)), rank)
        obj.label = label
        obj.rank = rank
        return obj

    def _latex(self, printer):
        return f"{self.label}_{self.rank}"

    def to_second_quant(self):
        creators = []
        annihilators = []
        holes = []
        parts = []
        sum_limits = []
        
        # Standard CC Textbook Alphabet
        occ_labels = ["i", "j", "k", "l", "m", "n"]
        vir_labels = ["a", "b", "c", "d", "e", "f"]
        
        for k in range(self.rank):
            # Pass ONLY the clean letter to the Dummy constructor!
            # SymPy tracks the uniqueness internally.
            i = Dummy(occ_labels[k], below_fermi=True)
            a = Dummy(vir_labels[k], above_fermi=True)
            
            creators.append(Fd(a))
            annihilators.append(F(i))
            
            holes.append(i)
            parts.append(a)
            
            sum_limits.append((i, 1, oo))
            sum_limits.append((a, 1, oo))
            
        op_string = Mul(*(creators + annihilators))
        amplitude = ClusterAmplitude(self.label.lower(), holes, parts)
        weight = Rational(1, factorial(self.rank)**2)
        
        return Sum(weight * amplitude * op_string, *sum_limits)


class FluctuationPotential(Expr):
    """
    Abstract representation of the normal-ordered Fluctuation Potential (\Phi or W_N).
    Generates fresh Dummy indices for each specific projection space to prevent collisions.
    """
    is_commutative = False

    def __new__(cls, label="W_N"):
        obj = Expr.__new__(cls, Symbol(str(label)))
        obj.label = label
        return obj

    def _latex(self, printer):
        return f"\\Phi"

    def to_second_quant(self, term="all"):
        """
        Translates the abstract potential into explicit normal-ordered strings.
        Allows for term-wise extraction (e.g., term="vvoo") for isolated derivations.
        """
        def get_dummies():
            i, j, k, l = [Dummy(x, below_fermi=True) for x in "ijkl"]
            a, b, c, d = [Dummy(x, above_fermi=True) for x in "abcd"]
            return (i, j, k, l), (a, b, c, d)

        terms = {}
        
        (i, j, _, _), _ = get_dummies()
        terms['scalar'] = Sum(Rational(-1, 2) * DoubleBarIntegral((i, j), (i, j)), (i, 1, oo), (j, -oo, oo))

        # 2. oooo
        (i, j, k, l), _ = get_dummies()
        terms['oooo'] = Sum(Rational(1, 4) * DoubleBarIntegral((i, j), (k, l)) * F(l) * F(k) * Fd(i) * Fd(j), (i, 1, oo), (j, 1, oo), (k, 1, oo), (l, 1, oo))

        # 3. oovv
        (i, j, _, _), (a, b, _, _) = get_dummies()
        terms['oovv'] = Sum(Rational(1, 4) * DoubleBarIntegral((i, j), (a, b)) * Fd(i) * Fd(j) * F(b) * F(a), (i, 1, oo), (j, 1, oo), (a, 1, oo), (b, 1, oo))

        # 4. vvoo
        (i, j, _, _), (a, b, _, _) = get_dummies()
        terms['vvoo'] = Sum(Rational(1, 4) * DoubleBarIntegral((a, b), (i, j)) * Fd(a) * Fd(b) * F(j) * F(i), (a, 1, oo), (b, 1, oo), (i, 1, oo), (j, 1, oo))

        # 5. ovov
        (i, j, _, _), (a, b, _, _) = get_dummies()
        terms['ovov'] = Sum(DoubleBarIntegral((a, i), (j, b)) * Fd(a) * F(j) * Fd(i) * F(b), (a, -oo, oo), (i, 1, oo), (j, 1, oo), (b, 1, oo))

        # 6. ovvv
        (i, _, _, _), (a, b, c, _) = get_dummies()
        terms['ovvv'] = Sum(Rational(1, 2) * DoubleBarIntegral((i, a), (b, c)) * Fd(i) * Fd(a) * F(c) * F(b), (i, 1, oo), (a, 1, oo), (b, 1, oo), (c, 1, oo))

        # 7. vvvo
        (i, _, _, _), (a, b, c, _) = get_dummies()
        terms['vvvo'] = Sum(Rational(1, 2) * DoubleBarIntegral((a, b), (i, c)) * Fd(a) * Fd(b) * F(c) * F(i), (a, 1, oo), (b, 1, oo), (i, 1, oo), (c, 1, oo))

        # 8. ooov
        (i, j, k, _), (a, _, _, _) = get_dummies()
        terms['ooov'] = Sum(Rational(1, 2) * DoubleBarIntegral((i, j), (k, a)) * F(k) * Fd(i) * F(a) * Fd(j), (i, 1, oo), (j, 1, oo), (k, 1, oo), (a, 1, oo))

        # 9. vooo
        (i, j, k, _), (a, _, _, _) = get_dummies()
        terms['vooo'] = Sum(Rational(1, 2) * DoubleBarIntegral((k, a), (i, j)) * F(j) * Fd(a) * F(i) * Fd(k), (k, 1, oo), (a, 1, oo), (i, 1, oo), (j, 1, oo))

        # 10. vvvv
        _, (a, b, c, d) = get_dummies()
        terms['vvvv'] = Sum(Rational(1, 4) * DoubleBarIntegral((a, b), (c, d)) * Fd(a) * Fd(b) * F(d) * F(c), (a, 1, oo), (b, 1, oo), (c, 1, oo), (d, 1, oo))

        if term == "all":
            return Add(*terms.values())
        elif term in terms:
            return terms[term]
        else:
            valid_keys = ", ".join(terms.keys())
            raise ValueError(f"Term '{term}' not recognized. Valid options: {valid_keys}")

def refresh_dummies(expr):
    """
    Finds all Sum objects in an expression and replaces their bound variables 
    with completely fresh Dummy clones, preserving all original assumptions.
    """
    if not expr.has(Sum):
        return expr

    replacements = {}
    
    # 1. Find every Sum object in the expression tree
    for summation in expr.find(Sum):
        # summation.variables gives you the indices being summed over (e.g., i, a)
        for var in summation.variables:
            
            # 2. Create a fresh clone using the same name and assumptions
            # .assumptions0 holds the above_fermi/below_fermi flags!
            new_var = Dummy(var.name, **var.assumptions0)
            
            # 3. Map the old variable to the new variable
            replacements[var] = new_var

    # 4. Use SymPy's xreplace to safely swap them out simultaneously
    return expr.xreplace(replacements)
