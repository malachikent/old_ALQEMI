from sympy.physics.secondquant import Commutator, F, Fd
from sympy.physics.quantum import AntiCommutator
from sympy import Sum, Add, Mul, init_printing, pprint
from alqemi.core.tools import  FermionAntiCommutator

init_printing()
def alqemi_reduce(expr):
    """
    Custom rule to intercept a Commutator and expand it 
    based on the number of operators (parity) in the strings.
    """
    # --- 1. BASE CASES & SYMPY HIERARCHY ---
    if not hasattr(expr, 'args') or not expr.args:
        return expr

    if not isinstance(expr, (AntiCommutator, Commutator, Sum, Add, Mul)):
        return expr

    if isinstance(expr, Sum):
        return Sum(alqemi_reduce(expr.args[0]), *expr.args[1:])

    if isinstance(expr, Add):
        return Add(*[alqemi_reduce(arg) for arg in expr.args])

    if isinstance(expr, Mul):
        return Mul(*[alqemi_reduce(arg) for arg in expr.args])

    if isinstance(expr, Commutator):
        left, right = expr.args

        # --- 1. NESTED COMMUTATOR RESOLUTION ---
        if left.has(Commutator): left = alqemi_reduce(left)
        if right.has(Commutator): right = alqemi_reduce(right)

        # --- 2. SUMMATION EXTRACTION (NEW!) ---
        # Pull summation limits completely outside the commutator
        if isinstance(left, Sum):
            return Sum(alqemi_reduce(Commutator(left.args[0], right)), *left.args[1:])
        if isinstance(right, Sum):
            return Sum(alqemi_reduce(Commutator(left, right.args[0])), *right.args[1:])

        # --- 3. ALGEBRAIC DISTRIBUTION ---
        if isinstance(left, Add):
            return Add(*[alqemi_reduce(Commutator(arg, right)) for arg in left.args])
        if isinstance(right, Add):
            return Add(*[alqemi_reduce(Commutator(left, arg)) for arg in right.args])

        # --- 4. SCALAR EXTRACTION (NEW!) ---
        # Pull commutative tensors (Amplitudes, Integrals, Fractions) outside!
        if isinstance(left, Mul):
            scalars = [arg for arg in left.args if arg.is_commutative]
            ops = [arg for arg in left.args if not arg.is_commutative]
            if scalars:
                return Mul(*scalars) * alqemi_reduce(Commutator(Mul(*ops) if ops else 1, right))
                
        if isinstance(right, Mul):
            scalars = [arg for arg in right.args if arg.is_commutative]
            ops = [arg for arg in right.args if not arg.is_commutative]
            if scalars:
                return Mul(*scalars) * alqemi_reduce(Commutator(left, Mul(*ops) if ops else 1))

        # --- 5. PREFACTOR SEPARATION & PARITY COUNTING ---
        l_args = list(left.args) if isinstance(left, Mul) else [left]
        r_args = list(right.args) if isinstance(right, Mul) else [right]

        l_ops = [arg for arg in l_args if isinstance(arg, (F, Fd))]
        r_ops = [arg for arg in r_args if isinstance(arg, (F, Fd))]
        
        l_prefactors = [arg for arg in l_args if not isinstance(arg, (F, Fd))]
        r_prefactors = [arg for arg in r_args if not isinstance(arg, (F, Fd))]

        left_count = len(l_ops)
        right_count = len(r_ops)
        
        left_is_even = (left_count % 2 == 0)
        right_is_even = (right_count % 2 == 0)

        global_prefactor = Mul(*(l_prefactors + r_prefactors))
        new_expression = 0

        # --- 6. HELGAKER PARITY RULES ---
        # Rule 1: [even, single]
        if right_count == 1 and left_is_even:
            for kk in range(1, left_count + 1):
                before = Mul(*l_ops[:kk-1]) if len(l_ops[:kk-1]) > 0 else 1
                center = l_ops[kk-1]
                after = Mul(*l_ops[kk:]) if len(l_ops[kk:]) > 0 else 1
                
                anticom = FermionAntiCommutator(center, r_ops[0]).doit()
                term = ((-1)**kk) * before * anticom * after 
                new_expression += term

        # Rule 2: [odd, single]
        elif right_count == 1 and not left_is_even:
            before = Mul(*l_ops[:-1]) if len(l_ops[:-1]) > 0 else 1
            after = Mul(*l_ops[1:]) if len(l_ops[1:]) > 0 else 1
            
            term1 = Commutator(before, r_ops[0]) * l_ops[-1]
            term2 = l_ops[0] * Commutator(after, r_ops[0])
            
            new_expression += alqemi_reduce(term1) + alqemi_reduce(term2)

        # Rule 3: [even, string]
        elif left_is_even:
            for kk in range(1, right_count + 1):
                before = Mul(*r_ops[:kk-1]) if len(r_ops[:kk-1]) > 0 else 1
                center = r_ops[kk-1]
                after = Mul(*r_ops[kk:]) if len(r_ops[kk:]) > 0 else 1
                
                left_op_string = Mul(*l_ops)
                term = before * alqemi_reduce(Commutator(left_op_string, center)) * after
                new_expression += term

        # Rule 4: [odd, even string]
        elif right_is_even:
            for kk in range(1, left_count + 1):
                before = Mul(*l_ops[:kk-1]) if len(l_ops[:kk-1]) > 0 else 1
                center = l_ops[kk-1]
                after = Mul(*l_ops[kk:]) if len(l_ops[kk:]) > 0 else 1
                
                right_op_string = Mul(*r_ops)
                term = before * alqemi_reduce(Commutator(center, right_op_string)) * after
                new_expression += term
                
        else:
            raise ValueError("ODD/ODD commutators not supported in this reduction scheme.")

        return (global_prefactor * new_expression).expand()

    if isinstance(expr, AntiCommutator):

        left, right = expr.args

        if left.has(AntiCommutator): left = alqemi_reduce(left)
        if right.has(AntiCommutator): right = alqemi_reduce(right)

        # --- SUMMATION EXTRACTION ---
        if isinstance(left, Sum):
            return Sum(alqemi_reduce(AntiCommutator(left.args[0], right)), *left.args[1:])
        if isinstance(right, Sum):
            return Sum(alqemi_reduce(AntiCommutator(left, right.args[0])), *right.args[1:])

        # --- ALGEBRAIC DISTRIBUTION ---
        if isinstance(left, Add):
            return Add(*[alqemi_reduce(AntiCommutator(arg, right)) for arg in left.args])
        if isinstance(right, Add):
            return Add(*[alqemi_reduce(AntiCommutator(left, arg)) for arg in right.args])

        # --- SCALAR EXTRACTION ---
        if isinstance(left, Mul):
            scalars = [arg for arg in left.args if arg.is_commutative]
            ops = [arg for arg in left.args if not arg.is_commutative]
            if scalars:
                return Mul(*scalars) * alqemi_reduce(AntiCommutator(Mul(*ops) if ops else 1, right))
                
        if isinstance(right, Mul):
            scalars = [arg for arg in right.args if arg.is_commutative]
            ops = [arg for arg in right.args if not arg.is_commutative]
            if scalars:
                return Mul(*scalars) * alqemi_reduce(AntiCommutator(left, Mul(*ops) if ops else 1))

        l_args = list(left.args) if isinstance(left, Mul) else [left]
        r_args = list(right.args) if isinstance(right, Mul) else [right]

        l_ops = [arg for arg in l_args if isinstance(arg, (F, Fd))]
        r_ops = [arg for arg in r_args if isinstance(arg, (F, Fd))]
        
        l_prefactors = [arg for arg in l_args if not isinstance(arg, (F, Fd))]
        r_prefactors = [arg for arg in r_args if not isinstance(arg, (F, Fd))]

        left_count = len(l_ops)
        right_count = len(r_ops)
        
        left_is_even = (left_count % 2 == 0)
        right_is_even = (right_count % 2 == 0)

        global_prefactor = Mul(*(l_prefactors + r_prefactors))
        new_expression = 0

        # --- HELGAKER PARITY RULES (AntiCommutator) ---
        if right_count == 1 and not left_is_even:
            for kk in range(1, left_count + 1):
                before = Mul(*l_ops[:kk-1]) if len(l_ops[:kk-1]) > 0 else 1
                center = l_ops[kk-1]
                after = Mul(*l_ops[kk:]) if len(l_ops[kk:]) > 0 else 1
                
                anticom = FermionAntiCommutator(center, r_ops[0]).doit()
                term = ((-1)**(kk-1)) * before * anticom * after 
                new_expression += term

        elif left_count == 1 and not right_is_even:
            for kk in range(1, right_count + 1):
                before = Mul(*r_ops[:kk-1]) if len(r_ops[:kk-1]) > 0 else 1
                center = r_ops[kk-1]
                after = Mul(*r_ops[kk:]) if len(r_ops[kk:]) > 0 else 1
                
                anticom = FermionAntiCommutator(l_ops[0], center).doit()
                term = ((-1)**(kk-1)) * before * anticom * after 
                new_expression += term

        elif not right_is_even and not left_is_even:
            for kk in range(1, right_count + 1):
                before = Mul(*r_ops[:kk-1]) if len(r_ops[:kk-1]) > 0 else 1
                center = r_ops[kk-1]
                after = Mul(*r_ops[kk:]) if len(r_ops[kk:]) > 0 else 1
                
                left_op_string = Mul(*l_ops)
                term = ((-1)**(kk-1)) * before * alqemi_reduce(AntiCommutator(left_op_string, center)) * after
                new_expression += term

        elif left_is_even:
            raise ValueError("EVEN/STRING anticommutators not supported in this reduction scheme.")
        elif right_is_even:
            raise ValueError("ODD/EVEN STRING anticommutators not supported in this reduction scheme.")
        else:
            raise ValueError("Commutator parity mapping failed.")
            
        return (global_prefactor * new_expression).expand()

    return expr


