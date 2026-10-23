"""Numerical oracle for ALQEMI tests.

Builds exact matrix representations of fermion creation/annihilation
operators on a small Fock space (Jordan-Wigner construction) and lets tests
compare a symbolic result against brute-force matrix algebra for EVERY
assignment of orbital indices to modes.

Why this exists: it is an independent referee. It never uses Wick's theorem
and never calls ALQEMI's own algebra -- it only *reads* ALQEMI's data
objects (Index, Op, OpString). If ALQEMI's symbolic answer and the matrices
disagree for any index assignment, the symbolic answer is wrong.

Orbital spaces: modes 0..n_occ-1 are occupied, n_occ..n_occ+n_virt-1 are
virtual. An OCC index ranges over occupied modes, VIRT over virtual modes,
GENERAL over all modes.
"""
from __future__ import annotations

import itertools
from typing import Dict, Iterable, Iterator, List

import numpy as np

from alqemi.core.indices import Index, Space
from alqemi.core.operators import Kind, Op, OpString

Assignment = Dict[Index, int]


class FockSpace:
    """Matrix fermion operators on n_occ + n_virt spin-orbitals."""

    def __init__(self, n_occ: int = 2, n_virt: int = 2):
        if n_occ < 1 or n_virt < 1:
            raise ValueError("need at least one occupied and one virtual mode")
        self.n_occ = n_occ
        self.n_virt = n_virt
        self.n_modes = n_occ + n_virt
        self.dim = 2 ** self.n_modes

        ident = np.eye(2)
        parity = np.diag([1.0, -1.0])        # Jordan-Wigner sign string
        lower = np.array([[0.0, 1.0], [0.0, 0.0]])  # |1> -> |0>

        self._ann: List[np.ndarray] = []
        for p in range(self.n_modes):
            mats = [parity] * p + [lower] + [ident] * (self.n_modes - p - 1)
            out = np.array([[1.0]])
            for m in mats:
                out = np.kron(out, m)
            self._ann.append(out)
        self._cre = [a.T.copy() for a in self._ann]  # real matrices: dagger = transpose

    # ---- index handling -------------------------------------------------
    def modes(self, space: Space) -> range:
        if space is Space.OCC:
            return range(0, self.n_occ)
        if space is Space.VIRT:
            return range(self.n_occ, self.n_modes)
        return range(0, self.n_modes)

    def assignments(self, indices: Iterable[Index]) -> Iterator[Assignment]:
        """Every way of assigning the given indices to modes allowed by their space."""
        unique: List[Index] = list(dict.fromkeys(indices))  # dedupe, keep order
        ranges = [self.modes(ix.space) for ix in unique]
        for combo in itertools.product(*ranges):
            yield dict(zip(unique, combo))

    # ---- operators ------------------------------------------------------
    @property
    def identity(self) -> np.ndarray:
        return np.eye(self.dim)

    def op(self, op: Op, assignment: Assignment) -> np.ndarray:
        mode = assignment[op.index]
        return self._cre[mode] if op.kind is Kind.CRE else self._ann[mode]

    def string(self, ops: OpString, assignment: Assignment) -> np.ndarray:
        out = self.identity
        for op in ops:
            out = out @ self.op(op, assignment)
        return out


# ---- matrix helpers (plain definitions, not identities) -----------------
def commutator(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a @ b - b @ a


def anticommutator(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a @ b + b @ a


def delta(p: Index, q: Index, assignment: Assignment) -> float:
    return 1.0 if assignment[p] == assignment[q] else 0.0


def first_mismatch(fock: FockSpace, indices, lhs, rhs):
    """Return the first assignment where lhs(asg) != rhs(asg), or None.

    lhs and rhs are callables taking an Assignment and returning a matrix.
    Returning the failing assignment (rather than True/False) tells you
    exactly which index coincidence breaks, e.g. {P: 0, Q: 1, R: 0}.
    """
    for asg in fock.assignments(indices):
        if not np.allclose(lhs(asg), rhs(asg)):
            return {str(k): v for k, v in asg.items()}
    return None
