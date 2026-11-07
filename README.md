# Economic cost of food-borne diseases in West Africa: country-level cost-of-illness model

R code and input data for the manuscript "A Systematic Review of the Economic Cost of Food-Borne Diseases in West Africa: A Country-Level Cost-of-Illness Model".

The model estimates the annual cost of food-borne diarrhoeal disease in 15 West African countries (Benin, Burkina Faso, Cabo Verde, Cote d'Ivoire, The Gambia, Ghana, Guinea, Guinea-Bissau, Liberia, Mali, Niger, Nigeria, Senegal, Sierra Leone and Togo). Illnesses and deaths come from the 2021 national estimates of the WHO Foodborne Disease Burden Epidemiology Reference Group (2026 edition, Lake et al., 2026). Costs are in 2023 US dollars. Premature death is valued by the human capital approach and, as an alternative, by the value of a statistical life (VSL). Uncertainty is propagated with a Monte Carlo simulation of 10,000 draws (seed 20261107).

## Folders

| Folder | Contents |
|---|---|
| `r/` | `fbd_model.R`, the full model, Monte Carlo simulation, scenarios, hypothesis tests and reconciliation |
| `data/` | Input files and a short source note for each (World Bank indicators, ILO employment ratios, WHO national estimates, United States reference values) |
| `tables/` | Output tables written by the model (country results, scenarios, Monte Carlo summary, hypothesis tests, reconciliation) |
| `figures/` | Figure 1, the PRISMA 2020 flow diagram |
| `results/` | `run_output.txt`, the printed R output of a full run |

## How to run

From the top-level folder of the repository:

```
Rscript run_all.R
```

The script uses base R only (no packages to install). It rewrites the files in `tables/` and `results/run_output.txt`. A full run takes a few seconds.

## Key files

- `tables/country_results.csv`: cost components by country (USD million, share of GDP, USD per person).
- `tables/scenarios.csv`: one-way scenario analysis for the human capital and VSL totals.
- `tables/monte_carlo_summary.csv`: median and 95% simulation interval for each cost component.
- `tables/h2_grid.csv` and `tables/h3_results.csv`: results for the two hypothesis tests (H₀₁ and H₀₂ in the manuscript).
- `tables/reconciliation.csv`: model results set against the figures reported by the included studies.

## Data sources

World Bank World Development Indicators (2023); International Labour Organization modelled estimates (2022); WHO Foodborne Disease Estimates, 2026 edition (Lake et al., 2026); unit costs from Baral et al. (2020). Each file in `data/` has a source note with the exact series used.
