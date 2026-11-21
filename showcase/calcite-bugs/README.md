# Calcite model: original vs fixed

A showcase of bugs in `models/calcite/calcite.phr`. Nothing in `models/` or
`datagen.py` is changed.

```
original/calcite.phr            the template, placeholders filled in, otherwise unchanged
fixes.py                        the four fixes, as text edits on the original
compare.py                      builds the variants, runs them, draws the figures
variants/<name>/calcite.phr     each fix alone, and all fixes ("# FIX n:" marks the edits)
figures/fix1-min.png ...        one fix vs the original, side by side
figures/all_fixes.png           the original, each fix alone and all fixes, side by side
output/<name>.tsv               PHREEQC's selected output of every run
```

In every figure the rows share one y axis, so a flat line in one column
next to a curve in another is a real difference.

All variants get the same input: the smoke-test values of `datagen.py`,
which are the median of every parameter in `models/calcite/config.cfg` with
the volumes, CH4 and calcite amounts from `derive()` (T = 42.5 °C, porosity
0.2, 100 L water, 100 L gas, CH4 58.15 atm, 21.66 mol reactive calcite).

## Run

```bash
pip install numpy matplotlib
python compare.py --phreeqc /path/to/phreeqc --database /path/to/phreeqc.dat
```

Without the options, it uses `executable` and `database` from
`models/calcite/config.cfg`. The six models run in parallel, each in 3 to
10 s. The figures were made with PHREEQC built from
`phreeqc-3.9.0-17591.tar.gz`. Its banner prints `PHREEQC_3.8.9`, so that is
the version shown in the figure titles.

## The fixes

| Fix | Where | Bug | Effect of the fix alone |
|---|---|---|---|
| 1 `fix1-min` | SR_H2, SR_AC, SR_FOR, IR_AC, IR_FOR, ACETO_H2, FORMATE_H2, METH_H2, X_SRB, X_IRB | `MIN()` is not a PHREEQC BASIC function. PHREEQC reads it as an undeclared array and it evaluates to 0 (checked: `PRINT MIN(3,5)` prints 0). So `lim = 0`, and `IF (moles > lim) THEN moles = lim` sets the reaction to 0. | The big one. Sulfate reduction, acetogenesis and formate formation run. H2 goes from 4.86 to 1.00 mol at day 180, pH rises from 6.79 to 7.01, and X_SRB stays near its seed instead of decaying to 1.7e-9 mol. |
| 2 `fix2-end` | IR_H2, IR_AC, IR_FOR, METH_H2 | No `END` after `SAVE moles`, so the program runs on into `900 moles = 0 / 910 SAVE moles`. The last `SAVE` wins. | Only IR_H2 runs, because it is the only one of the four without `MIN`, and only for its 1e-8 mol (see issue 1 below). The others still hit fix 1's bug. |
| 3 `fix3-cushion` | X_SRB, X_IRB lines 260 and 270 | The 0.9 cushion is applied twice, so the growth limit is 0.81 instead of 0.9. | None: `MIN` still makes the limit 0. The bug only matters once fix 1 is in. |
| 4 `fix4-n2` | USER_PUNCH line 55 | `GAS("Ntg(g)")` is counted twice in `n_gas`. | Only the reported gas volume: 98.02 to 97.65 L. The chemistry is the same. |

**All fixes together:** H2 ends at 0.66 mol instead of 1.00 mol with fix 1
alone. Fix 2 lets METH_H2 run, which needs fix 1 as well, and adds IR_AC and
IR_FOR.

Not changed: SR_H2 has `200 SAVE moles` twice. It's harmless; the second
line replaces the first.

## Issues the fixes expose (not fixed)

Fixing these takes model decisions, not just bug fixes, so they are left as
they are. They make the fixed runs behave in ways that are not physical:

1. **`-m` caps the total progress of every rate.** PHREEQC stops a kinetic
   reaction when its reactant `m` reaches 0, and these rates are written as
   if they could run forever. ACETO_H2 (`-m 1`) uses its whole 1 mol in the
   first 0.6 days, which takes 4 mol of H2. FORMATE_H2 (`-m 1`) is used up
   by day 3. IR_H2, IR_AC and IR_FOR (`-m 1e-8`) are used up in the first
   step, so iron reduction can never exceed 1e-8 mol. This explains the H2
   drop from 4.9 to 0.1 mol in 3 days and the flat iron-reduction lines.
2. **FORMATE_H2 is not redox-balanced** for PHREEQC's element bookkeeping.
   Hdg and Formate are their own elements, so `HCO3- -1 Hdg -1 Formate- 1
   H2O 1` adds 1 H, removes 1 C and 2 O at no charge change: +1 valence per
   mol. When FORMATE_H2 is used up (day 3), pe jumps from -3.4 to +12.4,
   all sulfide is oxidized to sulfate, and Fe(2) precipitates as Fe(OH)3
   (5e-5 to 6.7e-4 mol). The reverse shuttle FORMATE_H2_EQ and ACETO_METH
   (-1 valence per mol) partly offset it later.
3. **IR_* report a catabolic rate for reactions that don't run.** `PUT(r_exec,
   201..203)` happens before the fall-through of fix 2's bug, and after the
   reactant is used up (issue 1). X_IRB sees a growth signal larger than its
   decay, wants to grow, is capped by the available acetate (about 0) and
   stays at exactly 1e-8 mol in every run.
