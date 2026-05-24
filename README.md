# digester-cfd-adm1

A Python project, I built to model how well a biogas digester 
actually works -compared to how well it *should* work.

The main idea: real digesters have "dead zones" (areas where mixing is poor), 
which means microbes don't get enough substrate and produce less methane than expected. 
This project uses ADM1 (a standard biogas model) combined with CFD mixing data 
to show how big that difference actually is.

---

## What's inside

- **adm1_baseline.py** — compares ideal vs real methane output based on dead zone fraction
- **sensitivity_analysis.py** — shows how temperature, HRT and feeding rate affect CH4
- **cfd_dead_zones.py** — models 5 mixing scenarios from perfect to severely poor
- **adm1_full.py** — simulates the full biological process over 60 days using ODEs

## How to run

```bash
pip install numpy scipy pandas matplotlib
cd python_model
python3 adm1_baseline.py
```

## Results

Outputs (plots + CSV files) are saved in the `results/` folder.

---

Built with Python 3.14 
