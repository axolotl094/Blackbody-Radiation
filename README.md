# Blackbody Radiation: Lab Report (UP 201, Experiment 5)

My write-up and analysis for the blackbody radiation experiment, done in the lab on 26 October 2026 with my group (Anurag, Krishna, Aarush and Ananth).

We used a tungsten bulb as the blackbody, a prism spectrophotometer (PASCO OS-8537) to scan the light, and the bulb's resistance to work out its temperature. Then we checked Stefan's law and Wien's displacement law.

## What's in here

```
.
├── blackbody_report.tex        # the report (LaTeX source)
├── blackbody_report.pdf        # compiled report, 8 pages
├── blackbody_analysis.py       # analysis + all the graphs
├── results.json                # every number quoted in the report
├── figures/                    # PNGs used by the report
│   ├── fig1_resistance_room_temperature.png
│   ├── fig2_blackbody_spectrum_theory.png
│   ├── fig3_intensity_vs_angle.png
│   ├── fig4_intensity_vs_wavelength.png
│   ├── fig5_stefan_law.png
│   └── fig6_wien_law.png
```

## What I found

| Quantity | Result |
|---|---|
| Cold filament resistance R₀ | 3.62 ± 0.09 Ω |
| Filament temperature (6, 7, 7.5, 8 V) | 1515, 1611, 1663, 1713 K |
| λmax (same voltages) | 677, 660, 651, 634 nm |
| Stefan's law slope (expected 4) | 4.43 ± 0.34 |
| Wien's constant b | 1.06 × 10⁻³ m·K (accepted 2.898 × 10⁻³) |

Stefan's law works out well. Wien's law only works qualitatively: the peak does move to shorter wavelength as the bulb gets hotter, and λmax·T stays roughly constant, but the value of b comes out about a third of the real one. My measured peaks are around 650 nm, whereas a ~1600 K body should peak near 1.8 µm. The most likely reason is that the prism and light sensor don't respond equally across wavelengths, and I didn't correct for that. The report discusses this in detail.

## Things to know before trusting the numbers

- **Angle units.** The data files label the angle column "rad", but the values only make sense as degrees, so I treated them as degrees.
- **Angle zero.** I assumed the scan started at the undeviated beam (θ₀ = 0). If that's wrong, the wavelengths shift (about 90 nm per degree). You can change `THETA0` and `ANGLE_SCALE` at the top of the script.
- **Room temperature** wasn't recorded, so I assumed T₀ = 300 K.
- **Only four voltages** (6, 7, 7.5, 8 V) were usable. The 6.5 V and 9 V readings were struck out in the notebook.
- The lab handout has a sign typo in R = R₀[1 − α(T − T₀)]. I used a plus sign.

## How to run it

1. Put the data files in the same folder as the script :
   `6volt.txt`, `7volt.txt`, `7_5volt.txt`, `8volt.txt`, `anglecallibration.txt`
2. Install the dependencies and run:
   ```bash
   pip install numpy scipy matplotlib
   python3 blackbody_analysis.py
   ```
   This regenerates everything in `figures/` and `results.json`.
3. Build the report:
   ```bash
   pdflatex blackbody_report.tex
   pdflatex blackbody_report.tex   # second pass for references
   ```

## A note on AI assistance

I used an AI assistant (Claude, by Anthropic) to help with this report. It helped me write the analysis script, produce the graphs, and draft the LaTeX and the theory sections. The experiment, the measurements and the lab notebook are ours, and I checked the results and the assumptions listed above myself. If you reuse any of this, please do your own checking too, and follow your course's policy on AI use.

## License / use

This is coursework, shared for reference. Please don't submit it as your own work.
