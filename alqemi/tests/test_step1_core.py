"""Step 1 tests: the symbol layer (Space, Index, Op, OpString) and the oracle.

These tests DEFINE the Step 1 API. They fail until Claude Code implements:

  alqemi/core/indices.py
    class Space(Enum): GENERAL, OCC, VIRT
        .overlaps(other) -> bool    # can two indices in these spaces coincide?
    @dataclass(frozen=True) class Index:
        name: str                   # non-empty, else ValueError
        space: Space = Space.GENERAL   # must be a Space, else TypeError
        __str__ -> name

  alqemi/core/operators.py
    class Kind(Enum): CRE, ANN
    @dataclass(frozen=True) class Op:
        kind: Kind
        index: Index
        .adjoint() -> Op            # CRE <-> ANN, same index
        __mul__(Op | OpString) -> OpString
        __str__ -> "a†_p" (CRE) or "a_p" (ANN)
    def cre(index) -> Op ; def ann(index) -> Op
    @dataclass(frozen=True) class OpString:
        ops: tuple[Op, ...] = ()    # list accepted, stored as tuple; non-Op -> TypeError
        __len__, __iter__
        __getitem__(int) -> Op ; __getitem__(slice) -> OpString
        __mul__(Op | OpString) -> OpString      # concatenation, order preserved
        .is_even -> bool            # len % 2 == 0 (empty string is even)
        .split(k) -> (OpString, OpString)       # ops[:k], ops[k:]; 0<=k<=len else ValueError
        .indices -> tuple[Index, ...]           # unique, in order of first appearance
        __str__ -> ops joined by a single space; empty -> "1"

No algebra (brackets, deltas) is tested here. The oracle section checks the
oracle against the fermion anticommutation relations and records every
MATHEMATICAL_CONTRACT.md example as a golden result for later steps.
"""
import dataclasses

import numpy as np
import pytest

from alqemi.core.indices import Index, Space
from alqemi.core.operators import Kind, Op, OpString, ann, cre
from alqemi.tests.oracle import (
    FockSpace,
    anticommutator,
    commutator,
    delta,
    first_mismatch,
)

# General indices used throughout (named as in the contract)
P, Q, R, S, M, N = (Index(x) for x in "PQRSMN")
# Space-restricted indices
i, j = Index("i", Space.OCC), Index("j", Space.OCC)
a, b = Index("a", Space.VIRT), Index("b", Space.VIRT)


def E(p, q):
    """Unit-operator string a†_p a_q (data only, no algebra)."""
    return cre(p) * ann(q)


# =========================================================================
# Space
# =========================================================================
class TestSpace:
    def test_members(self):
        assert {s.name for s in Space} == {"GENERAL", "OCC", "VIRT"}

    @pytest.mark.parametrize(
        "s1, s2, expected",
        [
            (Space.GENERAL, Space.GENERAL, True),
            (Space.GENERAL, Space.OCC, True),
            (Space.GENERAL, Space.VIRT, True),
            (Space.OCC, Space.OCC, True),
            (Space.VIRT, Space.VIRT, True),
            (Space.OCC, Space.VIRT, False),  # this is what will make delta_ia = 0 later
        ],
    )
    def test_overlaps_is_symmetric(self, s1, s2, expected):
        assert s1.overlaps(s2) is expected
        assert s2.overlaps(s1) is expected


# =========================================================================
# Index
# =========================================================================
class TestIndex:
    def test_default_space_is_general(self):
        assert Index("p").space is Space.GENERAL

    def test_equality_is_by_value(self):
        assert Index("p") == Index("p", Space.GENERAL)
        assert Index("i", Space.OCC) == Index("i", Space.OCC)

    def test_same_name_different_space_are_different(self):
        assert Index("i", Space.OCC) != Index("i", Space.VIRT)
        assert Index("p") != Index("p", Space.OCC)

    def test_hashable_and_usable_in_sets(self):
        # Term.summed will be a frozenset[Index], so this must work.
        assert len({Index("p"), Index("p"), Index("q")}) == 2

    def test_immutable(self):
        p = Index("p")
        with pytest.raises(dataclasses.FrozenInstanceError):
            p.name = "q"

    def test_str_is_name(self):
        assert str(Index("i", Space.OCC)) == "i"

    def test_rejects_empty_name(self):
        with pytest.raises(ValueError):
            Index("")

    def test_rejects_non_space(self):
        with pytest.raises(TypeError):
            Index("p", "occ")


# =========================================================================
# Op
# =========================================================================
class TestOp:
    def test_helpers(self):
        assert cre(P) == Op(Kind.CRE, P)
        assert ann(P) == Op(Kind.ANN, P)

    def test_adjoint_flips_kind_keeps_index(self):
        assert cre(P).adjoint() == ann(P)
        assert ann(i).adjoint() == cre(i)

    def test_adjoint_is_involution(self):
        assert cre(a).adjoint().adjoint() == cre(a)

    def test_equality_and_hash(self):
        assert cre(P) == cre(Index("P"))
        assert cre(P) != ann(P)
        assert cre(P) != cre(Q)
        assert len({cre(P), cre(P), ann(P)}) == 2

    def test_immutable(self):
        with pytest.raises(dataclasses.FrozenInstanceError):
            cre(P).index = Q

    def test_str(self):
        assert str(cre(P)) == "a†_P"
        assert str(ann(i)) == "a_i"

    def test_mul_builds_opstring_in_order(self):
        s = cre(P) * ann(Q)
        assert isinstance(s, OpString)
        assert s.ops == (cre(P), ann(Q))


# =========================================================================
# OpString
# =========================================================================
class TestOpString:
    def test_empty(self):
        s = OpString()
        assert len(s) == 0
        assert s.is_even
        assert str(s) == "1"

    def test_list_input_stored_as_tuple(self):
        s = OpString([cre(P), ann(Q)])
        assert isinstance(s.ops, tuple)

    def test_rejects_non_op(self):
        with pytest.raises(TypeError):
            OpString((cre(P), "a_Q"))

    def test_len_and_iteration_preserve_order(self):
        s = cre(P) * ann(Q) * cre(R)
        assert len(s) == 3
        assert list(s) == [cre(P), ann(Q), cre(R)]

    @pytest.mark.parametrize("n, even", [(1, False), (2, True), (3, False), (4, True)])
    def test_parity(self, n, even):
        ops = [cre(P), ann(Q), cre(R), ann(S)][:n]
        assert OpString(ops).is_even is even

    def test_getitem_int_and_slice(self):
        s = cre(P) * ann(Q) * cre(R) * ann(S)
        assert s[0] == cre(P)
        assert s[-1] == ann(S)
        assert s[1:3] == ann(Q) * cre(R)
        assert isinstance(s[1:3], OpString)

    def test_concatenation_does_not_commute(self):
        # a†_P a_Q is a different string from a_Q a†_P; no reordering happens here.
        assert cre(P) * ann(Q) != ann(Q) * cre(P)

    def test_concatenation_with_strings_and_ops(self):
        s = (cre(P) * ann(Q)) * (cre(R) * ann(S))
        assert s.ops == (cre(P), ann(Q), cre(R), ann(S))
        assert (s * cre(M)).ops[-1] == cre(M)

    def test_concatenation_is_associative(self):
        x, y, z = E(P, Q), E(R, S), E(M, N)
        assert (x * y) * z == x * (y * z)

    @pytest.mark.parametrize("k", [0, 1, 2, 3, 4])
    def test_split_recombines(self, k):
        s = cre(P) * ann(Q) * cre(R) * ann(S)
        left, right = s.split(k)
        assert len(left) == k
        assert left * right == s

    @pytest.mark.parametrize("k", [-1, 5])
    def test_split_out_of_range(self, k):
        with pytest.raises(ValueError):
            (cre(P) * ann(Q) * cre(R) * ann(S)).split(k)

    def test_indices_unique_in_first_appearance_order(self):
        s = cre(Q) * ann(P) * cre(P) * ann(R)
        assert s.indices == (Q, P, R)

    def test_hashable_and_immutable(self):
        s = E(P, Q)
        assert len({s, E(P, Q), E(Q, P)}) == 2
        with pytest.raises(dataclasses.FrozenInstanceError):
            s.ops = ()

    def test_str(self):
        assert str(cre(P) * ann(Q)) == "a†_P a_Q"


# =========================================================================
# Oracle self-checks: the referee must obey the fermion relations
# =========================================================================
@pytest.fixture(scope="module")
def fock():
    return FockSpace(n_occ=2, n_virt=2)


class TestOracleAxioms:
    def test_assignments_respect_spaces(self, fock):
        asgs = list(fock.assignments([i, a, P]))
        assert len(asgs) == 2 * 2 * 4
        assert all(asg[i] in fock.modes(Space.OCC) for asg in asgs)
        assert all(asg[a] in fock.modes(Space.VIRT) for asg in asgs)

    def test_repeated_index_assigned_once(self, fock):
        assert len(list(fock.assignments([P, P, Q]))) == 4 * 4

    def test_cre_ann_anticommutator_is_delta(self, fock):
        # [a_P, a†_Q]_+ = delta_PQ
        bad = first_mismatch(
            fock, [P, Q],
            lambda g: anticommutator(fock.op(ann(P), g), fock.op(cre(Q), g)),
            lambda g: delta(P, Q, g) * fock.identity,
        )
        assert bad is None, bad

    def test_cre_cre_anticommutator_is_zero(self, fock):
        bad = first_mismatch(
            fock, [P, Q],
            lambda g: anticommutator(fock.op(cre(P), g), fock.op(cre(Q), g)),
            lambda g: 0 * fock.identity,
        )
        assert bad is None, bad

    def test_ann_ann_anticommutator_is_zero(self, fock):
        bad = first_mismatch(
            fock, [P, Q],
            lambda g: anticommutator(fock.op(ann(P), g), fock.op(ann(Q), g)),
            lambda g: 0 * fock.identity,
        )
        assert bad is None, bad

    def test_occ_virt_anticommutator_is_zero(self, fock):
        # [a_i, a†_a]_+ = delta_ia = 0 because i and a can never coincide
        bad = first_mismatch(
            fock, [i, a],
            lambda g: anticommutator(fock.op(ann(i), g), fock.op(cre(a), g)),
            lambda g: 0 * fock.identity,
        )
        assert bad is None, bad

    def test_oracle_detects_a_wrong_answer(self, fock):
        # Guard against a referee that passes everything: sign-flipped
        # contract example 1 must be caught.
        bad = first_mismatch(
            fock, [P, Q, R],
            lambda g: commutator(fock.op(cre(P), g), fock.string(E(Q, R), g)),
            lambda g: +delta(P, R, g) * fock.op(cre(Q), g),
        )
        assert bad is not None


# =========================================================================
# Golden results: every MATHEMATICAL_CONTRACT.md example, checked by matrices.
# Later steps will assert ALQEMI's symbolic output against these same RHS.
# =========================================================================
class TestContractExamples:
    def test_example_1(self, fock):
        # [a†_P, a†_Q a_R] = -delta_PR a†_Q
        bad = first_mismatch(
            fock, [P, Q, R],
            lambda g: commutator(fock.op(cre(P), g), fock.string(E(Q, R), g)),
            lambda g: -delta(P, R, g) * fock.op(cre(Q), g),
        )
        assert bad is None, bad

    def test_example_2(self, fock):
        # [a_P, a†_Q a_R] = delta_PQ a_R
        bad = first_mismatch(
            fock, [P, Q, R],
            lambda g: commutator(fock.op(ann(P), g), fock.string(E(Q, R), g)),
            lambda g: delta(P, Q, g) * fock.op(ann(R), g),
        )
        assert bad is None, bad

    def test_example_3(self, fock):
        # [E_PQ, E_RS] = delta_QR E_PS - delta_PS E_RQ
        bad = first_mismatch(
            fock, [P, Q, R, S],
            lambda g: commutator(fock.string(E(P, Q), g), fock.string(E(R, S), g)),
            lambda g: delta(Q, R, g) * fock.string(E(P, S), g)
            - delta(P, S, g) * fock.string(E(R, Q), g),
        )
        assert bad is None, bad

    def test_example_4(self, fock):
        # [E_PQ, E_RS E_MN] = d_QR E_PS E_MN - d_PS E_RQ E_MN + d_QM E_RS E_PN - d_PN E_RS E_MQ
        st = fock.string
        bad = first_mismatch(
            fock, [P, Q, R, S, M, N],
            lambda g: commutator(st(E(P, Q), g), st(E(R, S) * E(M, N), g)),
            lambda g: delta(Q, R, g) * st(E(P, S) * E(M, N), g)
            - delta(P, S, g) * st(E(R, Q) * E(M, N), g)
            + delta(Q, M, g) * st(E(R, S) * E(P, N), g)
            - delta(P, N, g) * st(E(R, S) * E(M, Q), g),
        )
        assert bad is None, bad

    def test_example_5(self, fock):
        # [E_PQ, [E_RS, E_MN]] = d_SM d_QR E_PN - d_SM d_PN E_RQ - d_RN d_QM E_PS + d_RN d_PS E_MQ
        st = fock.string
        bad = first_mismatch(
            fock, [P, Q, R, S, M, N],
            lambda g: commutator(st(E(P, Q), g), commutator(st(E(R, S), g), st(E(M, N), g))),
            lambda g: delta(S, M, g) * delta(Q, R, g) * st(E(P, N), g)
            - delta(S, M, g) * delta(P, N, g) * st(E(R, Q), g)
            - delta(R, N, g) * delta(Q, M, g) * st(E(P, S), g)
            + delta(R, N, g) * delta(P, S, g) * st(E(M, Q), g),
        )
        assert bad is None, bad


# =========================================================================
# The contract's general identities, checked on random matrices (any A, B_k).
# These justify that the only rule "condition" is the parity of n.
# =========================================================================
class TestContractIdentities:
    @pytest.fixture(scope="class")
    def mats(self):
        rng = np.random.default_rng(2026)
        return [rng.standard_normal((8, 8)) for _ in range(6)]

    @staticmethod
    def _prod(ms):
        out = np.eye(ms[0].shape[0])
        for m in ms:
            out = out @ m
        return out

    @pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
    def test_eq_1_8_5_any_n(self, mats, n):
        A, B = mats[0], mats[1 : 1 + n]
        rhs = sum(self._prod(B[:k] + [commutator(A, B[k])] + B[k + 1 :]) for k in range(n))
        assert np.allclose(commutator(A, self._prod(B)), rhs)

    @pytest.mark.parametrize("n", [2, 4])
    def test_eq_1_8_7_even_n(self, mats, n):
        A, B = mats[0], mats[1 : 1 + n]
        rhs = sum((-1) ** k * self._prod(B[:k] + [anticommutator(A, B[k])] + B[k + 1 :]) for k in range(n))
        assert np.allclose(commutator(A, self._prod(B)), rhs)

    @pytest.mark.parametrize("n", [1, 3, 5])
    def test_eq_1_8_9_odd_n(self, mats, n):
        A, B = mats[0], mats[1 : 1 + n]
        rhs = sum((-1) ** k * self._prod(B[:k] + [anticommutator(A, B[k])] + B[k + 1 :]) for k in range(n))
        assert np.allclose(anticommutator(A, self._prod(B)), rhs)

    @pytest.mark.parametrize("n", [3])
    def test_parity_condition_is_necessary(self, mats, n):
        # The alternating-anticommutator sum gives the COMMUTATOR only for even n
        # and the ANTICOMMUTATOR only for odd n. Swapping them must fail; this is
        # why Rule.matches() checks parity before applying 1.8.7 / 1.8.9.
        A, B = mats[0], mats[1 : 1 + n]
        alt = sum((-1) ** k * self._prod(B[:k] + [anticommutator(A, B[k])] + B[k + 1 :]) for k in range(n))
        assert not np.allclose(commutator(A, self._prod(B)), alt)

    def test_eq_1_8_8_both_forms(self, mats):
        A, B1, B2 = mats[:3]
        lhs = anticommutator(A, B1 @ B2)
        assert np.allclose(lhs, commutator(A, B1) @ B2 + B1 @ anticommutator(A, B2))
        assert np.allclose(lhs, anticommutator(A, B1) @ B2 - B1 @ commutator(A, B2))
