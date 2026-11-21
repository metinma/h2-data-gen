# Calcite model: original vs fixed

A showcase of bugs in `models/calcite/calcite.phr`. Nothing in `models/` or
`datagen.py` is changed.

```
original/calcite.phr     the template, placeholders filled in, otherwise unchanged
fixed/calcite.phr        the same file with the fixes below, each marked "# FIX:"
compare.py               runs both and draws calcite_comparison.png
calcite_comparison.png   the figure: original left, fixed right, same y axis per row
output/*.tsv             PHREEQC's selected output of both runs
```

Both files get the same input: the smoke-test values of `datagen.py`, which
are the median of every parameter in `models/calcite/config.cfg` with the
volumes, CH4 and calcite amounts from `derive()` (T = 42.5 °C, porosity 0.2,
100 L water, 100 L gas, CH4 58.15 atm, 21.66 mol reactive calcite).

## Run

```bash
pip install numpy matplotlib
python compare.py --phreeqc /path/to/phreeqc --database /path/to/phreeqc.dat
```

Without the options, it uses `executable` and `database` from
`models/calcite/config.cfg`. Each model takes about 3 to 10 s. The figure
was made with PHREEQC built from `phreeqc-3.9.0-17591.tar.gz`. Its banner
prints `PHREEQC_3.8.9`, so that is the version shown in the figure title.

## Bugs fixed

| # | Where | Bug | Effect |
|---|---|---|---|
| 1 | SR_H2, SR_AC, SR_FOR, IR_AC, IR_FOR, ACETO_H2, FORMATE_H2, METH_H2, X_SRB, X_IRB | `MIN()` is not a PHREEQC BASIC function. PHREEQC reads it as an undeclared array and it evaluates to 0 (checked: `PRINT MIN(3,5)` prints 0). | `lim = 0`, so `IF (moles > lim) THEN moles = lim` sets every one of these reactions to 0, and biomass can never grow. |
| 2 | IR_H2, IR_AC, IR_FOR, METH_H2 | No `END` after `SAVE moles`, so the program runs on into `900 moles = 0 / 910 SAVE moles`. The last `SAVE` wins. | These reactions are always 0, even where `lim` doesn't use `MIN`. |
| 3 | X_SRB, X_IRB lines 260 and 270 | The 0.9 cushion is applied twice. | The growth limit is 0.81 instead of 0.9 (hidden by bug 1). |
| 4 | USER_PUNCH line 55 | `GAS("Ntg(g)")` is counted twice in `n_gas`. | `Gas_VOL` is too high by the N2 volume. |

Not changed: SR_H2 has `200 SAVE moles` twice. It's harmless; the second
line replaces the first.

What the figure shows: in the original, sulfate reduction, iron reduction,
acetogenesis, formate formation and methanogenesis are all exactly 0. H2
barely moves (4.95 to 4.86 mol in 180 days). In the fixed model they run.

## Issues the fixes expose (not fixed)

Fixing these takes model decisions, not just bug fixes, so they are left as
they are. They make the fixed run behave in ways that are not physical:

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
   201..203)` happens before the fall-through of bug 2 in the original, and
   after the reactant is used up (issue 1) in the fixed model. X_IRB sees a
   growth signal larger than its decay, wants to grow, is capped by the
   available acetate (about 0) and stays at exactly 1e-8 mol in both runs.
