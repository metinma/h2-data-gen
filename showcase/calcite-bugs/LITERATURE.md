# Calcite model vs the Mura and Muller papers

Which experiment `models/calcite/calcite.phr` is built on, and how far it is
from it.

Papers checked:

- **Mura 2024a**: Mura et al. (2024), *Experimental simulation of H2
  coinjection via a high-pressure reactor with natural gas in a low-salinity
  deep aquifer used for current underground gas storage*, Front. Microbiol.
  15:1439866. Site **Ab_L**.
- **Mura 2024b**: Mura et al. (2024), *Comparative study of three H2
  geological storages in deep aquifers simulated in high-pressure reactors*,
  Int. J. Hydrogen Energy 63, 330–345. Sites Pb_T_1, **Pb_J_11**, Pb_C_5.
- **Mura 2025**: Mura et al. (2025), *Microbiological and Geochemical
  Perspectives on Sustainable Dihydrogen Storage in Deep Aquifers*, Environ.
  Sci. Technol. Site Pb_T_5.
- **Muller 2024**: Muller et al. (2024), *Evidencing the influence of
  temperature and mineralogy on microbial competition for hydrogen
  consumption*, Int. J. Hydrogen Energy 82, 1101–1113.

## Short answer

The template is built on the **Ab_L experiment of Mura 2024a**: 36 °C,
60 bar, 2% H2 in natural gas, the same four H2 reactions with the same
stoichiometry (including the 0.93 / 0.07 H2S / HS- split, which only that
paper uses), and the same dissolved iron. But it is not a model of that
experiment:

- The **water** is site **Pb_J_11** of Mura 2024b, copied exactly. Ab_L's
  water is about 10 to 100 times more dilute.
- The **rock** is 100% calcite. Ab_L is 80% quartz and 12% calcite. No site
  in these papers has a calcite rock. The template's title, "Calcite Midyat
  formation", names a formation in Turkey that is in none of the papers.
- The **CO2** is about 6% of the gas, plus 1.6% N2. Every Mura experiment
  used 99% CH4 and 1% CO2.
- Half the reactions (iron reduction, biomass, SR on acetate and formate,
  acetoclastic methanogenesis, the formate shuttle) and **all the rate
  constants** are in none of the four papers. Mura 2024a and 2024b compute
  extents of reaction from the measured data; they fit no rate constants.
  Mura 2025 fits a kinetic model, but for another site (53 °C, 107.6 bar,
  Thermoddem database), and its parameters are in its Supporting
  Information, which I don't have.
- The **results** don't match Ab_L, with or without the bug fixes (below).

Muller 2024 uses calcite powder, but at 2 bar with an 80/20 H2/CO2 gas and a
Gompertz fit of H2 consumption. Nothing in the template comes from it.

## Inputs

| | Template | Ab_L (Mura 2024a) | Pb_J_11 (Mura 2024b) |
|---|---|---|---|
| Temperature | 36 °C (comment; config samples 25–60) | **36 °C** | 35–36 °C |
| Pressure | `-pressure 60` (atm; the comment says 60 bar = 59.2 atm) | **60 bar target, 62 bar measured** | 80–85 bar |
| H2 | 1.1843 atm, "always 2%" | **2.2%** (one injection on day 9) | 8.9% (day 43) |
| CO2 | 3.68 atm (6.2%) | 1% | 1% |
| N2 | 0.92 atm (1.6%) | none | none |
| HCO3- | 34.45 mM | 3.44 mM | **34.45 mM** |
| Na | 62.78 mM | 0.50 mM | **62.782 mM** |
| Cl | 11.62 mM | 0.23 mM | **11.62 mM** |
| SO4 | 13.03 mM | 0.09 mM | **13.03 mM** |
| Ca | 1.24 mM | 1.37 mM | **1.24 mM** |
| Mg | 1.57 mM (config: 0.1 mM median) | 0.30 mM | **1.57 mM** |
| Fe(2) | 6.5e-6 mol/kg | **6.5e-6 mol/kg** (initial, measured) | 0.08 µM |
| Rock | 100% calcite | 80% quartz, 12% calcite | quartz and clays |
| Duration | 180 days | 105 days (96 after H2) | 78 days |

Bold marks a value equal to the template's.

## Reactions

| Template | Mura 2024a (Ab_L) | Mura 2024b | Mura 2025 |
|---|---|---|---|
| SR_H2: SO4 + 4 H2 + 1.93 CO2 → 0.93 H2S + 0.07 HS- + 1.93 HCO3- + 2.07 H2O | **Eq. 1, identical** | Eq. 2 with ν = 0.85 / 0.67 / 0.56 | SR on H2 and on formate |
| ACETO_H2: CO2 + HCO3- + 4 H2 → acetate + 3 H2O | **Eq. 2, identical** | Eq. 3 (printed with 2 H2O, which doesn't balance) | yes |
| FORMATE_H2: HCO3- + H2 → formate + H2O | **Eq. 3** | Eq. 4 | yes |
| METH_H2: CO2 + 4 H2 → CH4 + 2 H2O | **Eq. 4** | Eq. 1 | on H2 and on formate |
| SR_AC, SR_FOR, ACETO_METH, FORMATE_H2_EQ | no | no | SR on formate only |
| IR_H2, IR_AC, IR_FOR, Fe(OH)3(a) 5e-5 mol | no | no | no (pyrrhotite and carbonate precipitation instead) |
| X_SRB, X_IRB biomass | no | no | methanogens and H2-SRB, GEDYM yields, variable death rate |
| Monod rates, "p_fit from MATLAB" | none (extent of reaction) | none (extent of reaction) | double Monod × pH factor × thermodynamic factor; values in the SI |
| Database | phreeqc.dat with Hdg and Mtg decoupled | phreeqc.dat for gas–liquid; Thermoddem with the slow C, N, S redox couples removed | modified Thermoddem, slow C, N, S redox removed |

On the database: Mura removes the slow sulfur redox couples, so sulfate and
sulfide can't equilibrate. The template keeps S(6)/S(-2) coupled. That is
why, in the fixed run, the unbalanced FORMATE_H2 (README, issue 2) can
oxidize all the sulfide to sulfate.

## Results vs Ab_L

Ab_L numbers are from Mura 2024a §3.1, §3.2 and §3.6. Template numbers are
from `output/original.tsv` and `output/all-fixes.tsv` (180 days, 1 m³ REV).

| | Ab_L | Template, original | Template, all fixes |
|---|---|---|---|
| H2 lost | 28% (0.106 → 0.076 mol in 96 days) | 2% | 87% (4.95 → 0.66 mol) |
| Sulfate | used up 6 days after H2 injection | unchanged | unchanged (SR reacts about 1e-6 of 0.9 mol) |
| Main H2 sink | methanogenesis, 65% (SR 9%, acetogenesis 14%, formate 12%) | none: all four H2 reactions are 0 | acetogenesis, about 92% (4 mol in 0.6 days, then capped by `-m 1`); methanogenesis about 8% |
| Formate | up to 1.1e-3 mol/kg (day 21) | 9.1e-4 mol/kg | 1.8e-4 mol/kg |
| Acetate | up to 1.1e-3 mol/kg | ~1e-12 mol/kg | ~2e-12 mol/kg |
| pH | in situ, modelled from samples | 6.79 | 7.04 |

The original model consumes almost no H2 because of the bugs. The fixed
model consumes the wrong amount through the wrong pathway:

- Sulfate reduction is scaled by a biomass seed of 1e-8 mol, so it never
  starts.
- METH_H2 has K = 0.0002 d⁻¹ ("very low K to deprioritize H2 → CH4") and a
  day-6 gate, but methanogenesis was the main sink in Ab_L.
- Acetogenesis has K = 0.2 d⁻¹ without biomass and runs until its `-m`
  runs out.

## What it would take to fit Ab_L

These are model changes, not bug fixes, so they are not made here:

1. Use Ab_L's water (Mura 2024a Table 1, or the measured day-0 molalities
   in §3.2), 1% CO2 and no N2, and quartz 80% / calcite 12% for the rock.
2. Run 105 days with H2 injected on day 9, the water top-up on day 17, and
   the gas mix re-injection on day 7, rather than H2 from day 0 for 180
   days.
3. Fit the rate constants to the Ab_L data, or take the double-Monod form
   and parameters from Mura 2025 (needs its Supporting Information), and
   give methanogens a biomass term.
4. Fix the issues in the README: the `-m` caps, the redox balance of the
   decoupled species, and the S(6)/S(-2) coupling (remove it as Mura does).
5. Drop or justify the reactions with no basis in these papers: iron
   reduction and its biomass, SR on acetate, and acetoclastic
   methanogenesis.
