# ALQEMI v3 progress

## Done (2026-10-22)
- Step 1 symbol layer: `core/indices.py` (`Space`, `Index`) and `core/operators.py` (`Kind`, `Op`, `cre`, `ann`, `OpString`), as frozen dataclasses with no sympy (`8cfa1cb`).
- Step 1 tests now cite HJO equation numbers and include golden (1S.2.4), whose result is not normal ordered (`edf25b9`). Docs updated (`c047037`). Removed the empty `alqemi/operators/` package.
- Pushed `claude_alqemi` to origin.

## In progress
- Nothing is half-done. Step 2 has no spec or tests yet.

## Tests
- `pytest alqemi`: 81 pass (70 in step 1, 11 legacy), 0 failing.

## Open questions
- ElementaryBracket: should [a_p, a†_q]_+ stay as an unevaluated scalar until stage 3? See DESIGN_DECISIONS.md, OPEN section.
- Default route: split the left string with (1.8.4), as HJO does, or the right string with (1.8.5)?
- GitHub redirects `malachikent/alqemi` to `old_ALQEMI`. Is that the intended remote?

## Next step
Settle the two DESIGN_DECISIONS.md OPEN items. Then write the Step 2 spec and tests (`alqemi/tests/test_step2_*.py`) before writing any Step 2 code.
