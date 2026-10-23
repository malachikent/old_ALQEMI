from sympy import Dummy, Symbol
from sympy.physics.secondquant import F, Fd, Commutator, wicks, evaluate_deltas
from sympy.physics.quantum import AntiCommutator
from sympy import Mul
from algebra.rules import alqemi_reduce
from core.tools import FermionAntiCommutator, Braket, StrictCommutator
from core.cluster import ClusterOperator, FluctuationPotential, refresh_dummies
from sympy import init_printing, pprint, sympify
from sympy.printing.pretty.pretty import pretty_print

init_printing(use_unicode=True, wrap_line=False)

##################################################
########## Initialization for playground #########
##################################################

p, q, r= [Dummy(x) for x in "pqr"]
s, m, n= [Dummy(x) for x in "smn"]
# 1. Define indices with explicit spaces
# i, j are occupied (holes)
i = Dummy('i', below_fermi=True)
j,k = [Dummy(x, below_fermi=True) for x in 'jk']

# a, b are virtual (particles)
a = Dummy('a', above_fermi=True)
b,c = [Dummy(x, above_fermi=True) for x in 'bc']

i_hat = Symbol('i_hat', below_fermi=True)
a_hat = Symbol('a_hat', above_fermi=True)
j_hat = Symbol('j_hat', below_fermi=True)
b_hat = Symbol('b_hat', above_fermi=True)


Phi = FluctuationPotential()
W = Phi.to_second_quant('all')

T1 = ClusterOperator('T',1).to_second_quant()
T2 = ClusterOperator('T',2).to_second_quant()
T3 = ClusterOperator('T',3).to_second_quant()

bra = 1
bra_singles = Fd(i_hat) * F(a_hat)
bra_doubles = Fd(i_hat) * Fd(j_hat) * F(a_hat) * F(b_hat)

# 2. Build an operator string
A = Fd(p)*Fd(q)*F(s)*F(r)*Fd(a)*Fd(b)*F(j)*F(i)
# A = Fd(i)*Fd(j)*F(b)*F(a)*Fd(r)*Fd(s)*F(q)*F(p)
# A = Fd(a)*F(p)

##################################################
################### Playground ###################
##################################################
# print('W expansion:')
# pprint(W)
# print('T1 cluster operator')
# pprint(T1)
# print('T2 cluster operator')
# pprint(T2)
# print('T3 cluster operator')
# pprint(T3)
#
# inner_comm = StrictCommutator(W,T2)
# # pprint(inner_comm)
#
# A = alqemi_reduce(inner_comm)
#
# print('alqemi reduce test on [W,T2]')
# pprint(A)
#
# expectation = Braket(bra,A)
expectation = wicks(A)
#
print('Braket test')
pprint(expectation.doit())
