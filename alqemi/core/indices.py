"""Orbital indices and the orbital spaces they range over."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Space(Enum):
    """Orbital space an index ranges over."""

    GENERAL = "general"
    OCC = "occ"
    VIRT = "virt"

    def overlaps(self, other: Space) -> bool:
        """Whether an index in this space can coincide with one in ``other``."""
        if self is Space.GENERAL or other is Space.GENERAL:
            return True
        return self is other


@dataclass(frozen=True)
class Index:
    """A named orbital index; its space is explicit, never parsed from the name."""

    name: str
    space: Space = Space.GENERAL

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name:
            raise ValueError(f"Index name must be a non-empty string, got {self.name!r}")
        if not isinstance(self.space, Space):
            raise TypeError(f"Index space must be a Space, got {self.space!r}")

    def __str__(self) -> str:
        return self.name
