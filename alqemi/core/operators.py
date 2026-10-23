"""Fermion creation/annihilation operators and ordered operator strings.

Data only: no algebra happens here. Concatenation preserves order exactly.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterator, Tuple, Union, overload

from alqemi.core.indices import Index


class Kind(Enum):
    CRE = "cre"
    ANN = "ann"


@dataclass(frozen=True)
class Op:
    """A single creation (a†_p) or annihilation (a_p) operator."""

    kind: Kind
    index: Index

    def __post_init__(self) -> None:
        if not isinstance(self.kind, Kind):
            raise TypeError(f"Op kind must be a Kind, got {self.kind!r}")
        if not isinstance(self.index, Index):
            raise TypeError(f"Op index must be an Index, got {self.index!r}")

    def adjoint(self) -> Op:
        return Op(Kind.ANN if self.kind is Kind.CRE else Kind.CRE, self.index)

    def __mul__(self, other: Union[Op, OpString]) -> OpString:
        return OpString((self,)) * other

    def __str__(self) -> str:
        dagger = "†" if self.kind is Kind.CRE else ""
        return f"a{dagger}_{self.index}"


def cre(index: Index) -> Op:
    return Op(Kind.CRE, index)


def ann(index: Index) -> Op:
    return Op(Kind.ANN, index)


@dataclass(frozen=True)
class OpString:
    """An ordered product of operators; the empty string is the identity."""

    ops: Tuple[Op, ...] = ()

    def __post_init__(self) -> None:
        ops = tuple(self.ops)
        for op in ops:
            if not isinstance(op, Op):
                raise TypeError(f"OpString elements must be Op, got {op!r}")
        object.__setattr__(self, "ops", ops)

    def __len__(self) -> int:
        return len(self.ops)

    def __iter__(self) -> Iterator[Op]:
        return iter(self.ops)

    @overload
    def __getitem__(self, key: int) -> Op: ...

    @overload
    def __getitem__(self, key: slice) -> OpString: ...

    def __getitem__(self, key):
        if isinstance(key, slice):
            return OpString(self.ops[key])
        return self.ops[key]

    def __mul__(self, other: Union[Op, OpString]) -> OpString:
        if isinstance(other, Op):
            return OpString(self.ops + (other,))
        if isinstance(other, OpString):
            return OpString(self.ops + other.ops)
        return NotImplemented

    @property
    def is_even(self) -> bool:
        return len(self.ops) % 2 == 0

    def split(self, k: int) -> Tuple[OpString, OpString]:
        """Return (ops[:k], ops[k:]); requires 0 <= k <= len."""
        if not 0 <= k <= len(self.ops):
            raise ValueError(f"split position {k} out of range [0, {len(self.ops)}]")
        return OpString(self.ops[:k]), OpString(self.ops[k:])

    @property
    def indices(self) -> Tuple[Index, ...]:
        """Distinct indices in order of first appearance."""
        return tuple(dict.fromkeys(op.index for op in self.ops))

    def __str__(self) -> str:
        return " ".join(str(op) for op in self.ops) if self.ops else "1"
