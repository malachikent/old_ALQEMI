# ALQEMI v3 — Design Decisions

Status key: SETTLED = agreed, copy into CLAUDE.md. OPEN = needs a decision.

## SETTLED (2026-10-22)

1. **No Wick's theorem, anywhere.** Do not use `sympy.physics.secondquant.Commutator`
   (its constructor reorders arguments and calls `wicks()` for two single fermion
   operators) or `wicks`. SymPy may be used only for final printing.
2. **Own immutable core.** Frozen dataclasses: every rewrite returns new objects, so
   nothing simplifies behind our back and every step can be logged. Exact
   `fractions.Fraction` coefficients.
3. **Orbital space lives on the Index object**, never parsed from the name.
   `Index(name, space)` with `Space.GENERAL | OCC | VIRT`; equality is by value
   (name AND space).
4. **Dummies are a property of a Term** (`Term.summed: frozenset[Index]`), not of an
   Index. Products of terms always freshen dummies through an explicit `IndexPool`.
   No global state.
5. **Identities are Rule objects** with `name`, `matches()`, `apply()`. The only
   condition on HJO 1.8.5 / 1.8.7 / 1.8.9 is the parity of the string being split
   (verified numerically for arbitrary A). Unsupported cases raise, never `quit()`.
6. **Pipeline:** Linearize → IdentityExpand (until only ElementaryBrackets remain) →
   DeltaEvaluate → Simplify. Kronecker deltas exist only from stage 3 on.
7. **Testing against a numerical oracle** (`tests/oracle.py`: Jordan–Wigner matrices,
   all index assignments respecting spaces), not against expression-tree equality.
   The oracle never uses Wick or ALQEMI's algebra.

## Step 1 (in progress): symbol layer
- `alqemi/core/indices.py`: `Space` (with `overlaps`), `Index`.
- `alqemi/core/operators.py`: `Kind`, `Op`, `cre`, `ann`, `OpString`.
- Spec + tests: `alqemi/tests/test_step1_core.py` (68 checks; validated against a
  throwaway reference implementation, and confirmed to catch planted bugs).

## OPEN
- **ElementaryBracket as an unevaluated scalar.** Proposal: elementary anticommutators
  [a_p, a†_q]_+ stay as unevaluated leaves through stages 1–2; linearity may move
  them as scalars (needed for nested brackets such as [E_PQ, [E_RS, E_MN]]), but no
  δ, δ_pp = 1, δ_ia = 0, or summation collapse before stage 3. Awaiting confirmation
  that this matches evaluation-order rule (3).
- **[odd, odd] commutators and anticommutators with an even string** do not reduce to
  deltas through HJO 1.8.4–1.8.9. Plan: raise `UnsupportedBracket`. To confirm against
  the Helgaker excerpts (not yet uploaded to the project).
