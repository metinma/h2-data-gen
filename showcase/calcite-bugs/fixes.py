"""The bug fixes for the calcite model, each as text edits on original/calcite.phr.

Each fix is a list of edits. An edit is a list of alternatives
(old, new, count): the first alternative whose `old` occurs exactly `count`
times is applied. Alternatives let a fix apply both on its own and on top of
another fix that rewrote the same lines (fix 3 after fix 1).
"""

MIN_NOTE = "# FIX 1: MIN() is not a PHREEQC BASIC function; it evaluates to 0, so lim = 0 and the reaction never runs"
END_NOTE = "# FIX 2: END was missing, so the program ran on into line 900 and always saved 0"


def _min(line_no, a, b, indent=""):
    """lim = 0.9 * MIN(a, b), written without MIN."""
    return (
        f"{indent}{MIN_NOTE}\n"
        f"{indent}{line_no} lim = {a}\n"
        f"{indent}{line_no + 1} IF ({b} < lim) THEN lim = {b}\n"
        f"{indent}{line_no + 2} lim = 0.9 * lim"
    )


FIXES = {
    "fix1-min": {
        "title": "Fix 1: MIN() replaced",
        "summary": "MIN() is not a PHREEQC BASIC function and evaluates to 0",
        "edits": [
            # X_SRB and X_IRB. The 0.9 is still applied twice here (that is fix 3).
            [("    260 lim = 0.9 * MIN(lim_ac, lim_dic)\n    270 lim = 0.9 * MIN(lim, lim_nh4)",
              f"    {MIN_NOTE}\n"
              "    260 lim = lim_ac\n"
              "    261 IF (lim_dic < lim) THEN lim = lim_dic\n"
              "    262 lim = 0.9 * lim\n"
              "    264 IF (lim_nh4 < lim) THEN lim = lim_nh4\n"
              "    270 lim = 0.9 * lim", 2)],
            [("170 lim = 0.9 * MIN(ac*vL/1.0, mFe/8.0)", _min(170, "ac*vL/1.0", "mFe/8.0"), 1)],       # IR_AC
            [("170 lim = 0.9 * MIN(fo*vL/1.0, mFe/2.0)", _min(170, "fo*vL/1.0", "mFe/2.0"), 1)],       # IR_FOR
            [("140 lim = 0.9*MIN(ac*vL, so4*vL)  # 0.9 leaves a small cushion for numerics",
              _min(140, "ac*vL", "so4*vL"), 1)],                                                       # SR_AC
            [("140 lim = 0.9*MIN((fo*vL)/4.0, (so4*vL))", _min(140, "(fo*vL)/4.0", "(so4*vL)"), 1)],  # SR_FOR
            [("150 lim = 0.9*MIN((h2*vL)/4.0, (so4*vL))", _min(150, "(h2*vL)/4.0", "(so4*vL)"), 1)],  # SR_H2
            [("150 lim = 0.9 * MIN((h2 * vL) / 4.0,(dic * vL))",
              _min(150, "(h2 * vL) / 4.0", "(dic * vL)"), 1)],                                         # ACETO_H2
            [("150 lim = 0.9*MIN((h2*vL)/1.0, (hco3*vL))", _min(150, "(h2*vL)/1.0", "(hco3*vL)"), 1)],  # FORMATE_H2
            [("160 lim = 0.9*MIN((h2*vL)/4.0, (co2*vL))", _min(160, "(h2*vL)/4.0", "(co2*vL)"), 1)],  # METH_H2
        ],
    },
    "fix2-end": {
        "title": "Fix 2: missing END added",
        "summary": "IR_H2, IR_AC, IR_FOR and METH_H2 run on into '900 moles = 0'",
        "edits": [
            [("200 SAVE moles\n900 moles = 0", f"200 SAVE moles\n{END_NOTE}\n210 END\n900 moles = 0", 3)],  # IR_*
            [("190 SAVE moles\n900 moles = 0", f"190 SAVE moles\n{END_NOTE}\n195 END\n900 moles = 0", 1)],  # METH_H2
        ],
    },
    "fix3-cushion": {
        "title": "Fix 3: 0.9 cushion applied once",
        "summary": "X_SRB and X_IRB apply the 0.9 cushion twice (0.81)",
        "edits": [
            [
                # On top of fix 1: drop the first cushion line.
                ("    262 lim = 0.9 * lim\n",
                 "    # FIX 3: the 0.9 cushion was applied twice (0.81); it is applied once now\n", 2),
                # On its own: drop the first 0.9 (MIN still makes lim 0 here).
                ("    260 lim = 0.9 * MIN(lim_ac, lim_dic)",
                 "    # FIX 3: the 0.9 cushion was applied twice (0.81); it is applied once now\n"
                 "    260 lim = MIN(lim_ac, lim_dic)", 2),
            ],
        ],
    },
    "fix4-n2": {
        "title": "Fix 4: N2 counted once",
        "summary": "USER_PUNCH counts Ntg(g) twice in the gas volume",
        "edits": [
            [('GAS("H2S(g)")+ GAS("Ntg(g)")\n', 'GAS("H2S(g)")   # FIX 4: Ntg(g) was counted twice\n', 1)],
        ],
    },
}


def apply(text, fix_names):
    for name in fix_names:
        for alternatives in FIXES[name]["edits"]:
            for old, new, count in alternatives:
                if text.count(old) == count:
                    text = text.replace(old, new)
                    break
            else:
                raise ValueError(f"{name}: no alternative matches: {alternatives[0][0]!r}")
    return text
