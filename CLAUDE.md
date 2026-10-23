## Core design (v3)
- No Wick's theorem. Never use sympy.physics.secondquant.Commutator or wicks; SymPy only for final printing.
- Own immutable core (frozen dataclasses); exact Fraction coefficients.
- Orbital space lives on Index(name, space) with Space.GENERAL|OCC|VIRT; never parsed from names.
- Dummies are Term.summed (frozenset[Index]); products freshen dummies via an explicit IndexPool. No globals.
- Each HJO identity is a Rule (name / matches / apply); the only condition on 1.8.5/1.8.7/1.8.9 is string parity.
  Unsupported cases raise, never quit().
- Pipeline: Linearize -> IdentityExpand -> DeltaEvaluate -> Simplify. Deltas exist only after stage 3.
- Tests check results against tests/oracle.py (matrix fermions), never expression-tree equality.
- Never edit tests/oracle.py or golden tests to make code pass.
