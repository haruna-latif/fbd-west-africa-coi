# =====================================================================
# Economic cost of food-borne diseases in 15 West African countries
# Country-level cost-of-illness model (human capital and VSL valuation)
# Burden year 2021 (WHO 2026 edition); monetary values in 2023 USD
# Base case: diarrhoeal illnesses and diarrhoeal deaths (same disease group)
# Run from the top-level folder with:  Rscript run_all.R
# =====================================================================
set.seed(20261107)
N      <- 10000      # Monte Carlo iterations
R_DISC <- 0.03       # discount rate for human capital mortality (Robinson et al. 2019 comparator)
RHO_C  <- 0.5        # assumed correlation of WHO uncertainty between countries
RHO_ID <- 0.5        # assumed correlation between the illness and death factors

d <- "data/"
wb  <- read.csv(paste0(d, "wb_2023.csv"))
lab <- read.csv(paste0(d, "labour_2023.csv"))
ppp <- read.csv(paste0(d, "gdp_ppp_2023.csv"))
inc <- read.csv(paste0(d, "income_cons_2021.csv"))
ill <- read.csv(paste0(d, "who_ill_2021.csv"))
dth <- read.csv(paste0(d, "who_dth_2021.csv"))
ig  <- read.csv(paste0(d, "who_ill_2021_groups.csv"))
dg  <- read.csv(paste0(d, "who_dth_2021_groups.csv"))
us  <- read.csv(paste0(d, "us_refs.csv"))
cpi <- function(y) us$value[us$item == "US_CPI_2010eq100" & us$year == y]
GNI_US <- us$value[us$item == "US_GNI_pc_PPP_current_intl" & us$year == 2023]
INFL_2015_2023 <- cpi(2023) / cpi(2015)            # US CPI, 2015 USD -> 2023 USD

# ---- country table -------------------------------------------------
cty <- merge(merge(merge(wb, lab[, c("iso3", "epr15")], by = "iso3"), ppp, by = "iso3"), inc, by = "iso3")
pick <- function(x, a, col) x[x$age == a, ][match(cty$iso3, x$iso3[x$age == a]), col]
for (k in c("low", "high", "mean")) {
  cty[[paste0("ill_", k)]] <- pick(ill, "all", k)
  cty[[paste0("dth_", k)]] <- pick(dth, "all", k)
}
cty$ill_u5 <- pick(ill, "u5", "mean");  cty$dth_u5 <- pick(dth, "u5", "mean")
cty$sh_ill_u5 <- cty$ill_u5 / cty$ill_mean
cty$sh_dth_u5 <- cty$dth_u5 / cty$dth_mean
gsum <- function(x, g, col) { x <- x[x$group %in% g, ]; tapply(x[[col]], x$iso3, sum)[cty$iso3] }
cty$sh_dia  <- gsum(ig, "DIA", "mean") / cty$ill_mean                 # diarrhoeal share of illnesses
cty$sh_inf  <- gsum(ig, c("DIA", "INF"), "mean") / cty$ill_mean       # diarrhoeal + invasive share
cty$dsh_dia <- gsum(dg, "DIA", "mean") / cty$dth_mean                 # diarrhoeal share of deaths
cty$ill_dia_low <- gsum(ig, "DIA", "low"); cty$ill_dia_high <- gsum(ig, "DIA", "high")
cty$dth_dia_low <- gsum(dg, "DIA", "low"); cty$dth_dia_high <- gsum(dg, "DIA", "high")
cty$epr <- cty$epr15 / 100
cty$emp_share <- cty$epr * (1 - cty$pop0_14_pct / 100)               # employed persons / population
cty$gdp <- cty$gdppc_usd * cty$pop
cty$pr  <- cty$gdppc_usd / cty$gdppc_ppp                             # market USD per international dollar
cty$cons <- ifelse(is.na(cty$cons_pct_2023), median(cty$cons_pct_2023, na.rm = TRUE), cty$cons_pct_2023) / 100
# VSL in 2023 international dollars (Robinson, Hammitt & O'Keeffe 2019)
vsl_eps <- function(e) pmax(20 * cty$gni_ppp, 160 * cty$gni_ppp * (cty$gni_ppp / GNI_US)^(e - 1))
cty$vsl_a <- vsl_eps(1.5)               # default: elasticity 1.5, floor 20 x GNI
cty$vsl_b <- 160 * cty$gni_ppp          # elasticity 1.0
cty$vsl_c <- 100 * cty$gni_ppp          # OECD-based ratio 100
cty$vsl_d <- vsl_eps(2.0)               # elasticity 2.0 (floor binds)
nc <- nrow(cty)

# ---- helpers -------------------------------------------------------
ann <- function(a, b, r) sum((1 + r)^-(a:b))           # PV of 1 per year, years a..b
# beta distribution with mean m whose 2.5th and 97.5th percentiles match lo and hi
fit_beta <- function(m, lo, hi) {
  f <- function(p) { a <- exp(p[1]); b <- exp(p[2])
    1e4 * (a / (a + b) - m)^2 + 1e2 * (qbeta(0.025, a, b) - lo)^2 + 1e2 * (qbeta(0.975, a, b) - hi)^2 }
  o <- optim(c(log(m * 20), log((1 - m) * 20)), f, control = list(reltol = 1e-12, maxit = 5000))
  c(a = exp(o$par[1]), b = exp(o$par[2])) }
# lognormal multiplier with mean 1 and SD chosen so that log(max/min) = 5.15 SD (assumption)
ln_sd <- function(lo, hi) (log(hi) - log(lo)) / 5.15

# ---- cost engine ---------------------------------------------------
# u = unit costs per episode in 2015 USD (scalars or one value per country)
cost <- function(ill_all, dth_all, p_care, p_hosp, u, d_out, d_in, wmult = 1, age_ad = 40,
                 r = R_DISC, caregiver = FALSE, inf = FALSE, dscope = "DIA",
                 gdppc = cty$gdppc_usd, net = FALSE, infl = INFL_2015_2023) {
  sh <- if (inf) cty$sh_inf else cty$sh_dia
  I  <- ill_all * sh
  n_h <- I * p_hosp; n_o <- I * (p_care - p_hosp)
  DMC  <- (n_o * u$dm_o + n_h * u$dm_i) * infl
  DNMC <- (n_o * u$dnm_o + n_h * u$dnm_i) * infl
  days <- p_hosp * d_in + (1 - p_hosp) * d_out           # calendar days of illness
  W <- gdppc / cty$emp_share * wmult; dayval <- W / 365  # value per calendar day
  Prod <- I * (1 - cty$sh_ill_u5) * cty$emp_share * days * dayval
  if (caregiver) Prod <- Prod + I * cty$sh_ill_u5 * cty$epr * days * dayval
  dths <- if (dscope == "DIA") dth_all * cty$dsh_dia else dth_all
  A_u5 <- ann(13, 62, r); A_ad <- ann(1, 64 - age_ad, r)
  d5 <- dths * cty$sh_dth_u5; da <- dths - d5
  MortHC <- (d5 * A_u5 + da * A_ad) * W * cty$epr * (if (net) 1 - cty$cons else 1)
  list(DMC = DMC, DNMC = DNMC, Prod = Prod, MortHC = MortHC,
       MortA = dths * cty$vsl_a * cty$pr, MortB = dths * cty$vsl_b * cty$pr,
       MortC = dths * cty$vsl_c * cty$pr, MortD = dths * cty$vsl_d * cty$pr,
       MortA_int = dths * cty$vsl_a)
}
tot <- function(x) { m <- x$DMC + x$DNMC + x$Prod
  c(DMC = sum(x$DMC), DNMC = sum(x$DNMC), Prod = sum(x$Prod), MortHC = sum(x$MortHC),
    TotalHC = sum(m + x$MortHC), VSL_A = sum(m + x$MortA), VSL_B = sum(m + x$MortB),
    VSL_C = sum(m + x$MortC), VSL_D = sum(m + x$MortD)) }

# ---- parameters ----------------------------------------------------
# Unit costs: Baral et al. 2020 Table 3, modelled weighted means by World Bank income group (2015 USD)
U_LIC  <- c(dm_o = 4.73, dnm_o = 1.63, dm_i = 24.40, dnm_i = 5.03)
U_LMIC <- c(dm_o = 6.17, dnm_o = 2.13, dm_i = 57.33, dnm_i = 11.81)
ulist <- function(lic, lmic) { v <- function(k) ifelse(cty$income == "LIC", lic[k], lmic[k])
  setNames(lapply(c("dm_o", "dnm_o", "dm_i", "dnm_i"), v), c("dm_o", "dnm_o", "dm_i", "dnm_i")) }
U <- ulist(U_LIC, U_LMIC)        # Cabo Verde (upper-middle income) takes the lower-middle-income values
# Baral Table 2 min-max across studies (one-way range and spread of the Monte Carlo multipliers)
RNG <- list(dm_o = c(0.92, 58.23), dnm_o = c(0.47, 18.87), dm_i = c(24.24, 341.27), dnm_i = c(1.18, 39.75))
SDU <- sapply(RNG, function(x) ln_sd(x[1], x[2]))
P <- list(p_care = 0.52, p_hosp = 0.03, d_out = 4, d_in = 4.1)
bc <- fit_beta(0.52, 0.33, 0.62); bh <- fit_beta(0.03, 0.013, 0.069)
cat("Beta check, care seeking 2.5/97.5%:", round(qbeta(c(.025, .975), bc["a"], bc["b"]), 3),
    " hospitalised:", round(qbeta(c(.025, .975), bh["a"], bh["b"]), 3), "\n")
SD_DAYS <- 0.20; SD_W <- 0.15
args_base <- c(list(ill_all = cty$ill_mean, dth_all = cty$dth_mean, u = U), P)
run <- function(...) do.call(cost, modifyList(args_base, list(...)))

# ---- 1. deterministic base case -----------------------------------
base <- run()
cat("\n1. BASE CASE, diarrhoeal illnesses and deaths (USD billion, 2023 prices; VSL in market USD)\n")
print(round(tot(base) / 1e9, 3))
cat("Illnesses, all hazards (million):", round(sum(cty$ill_mean) / 1e6, 2),
    " costed diarrhoeal (million):", round(sum(cty$ill_mean * cty$sh_dia) / 1e6, 2), "\n")
cat("Deaths, all hazards:", round(sum(cty$dth_mean)), " diarrhoeal deaths:", round(sum(cty$dth_mean * cty$dsh_dia)), "\n")
cat("VSL-A mortality in international dollars (bn):", round(sum(base$MortA_int) / 1e9, 2), "\n")
gdp <- sum(cty$gdp)
cat("Regional GDP (USD bn):", round(gdp / 1e9, 1), "\n")
cat("Share of GDP (%): HC total", round(100 * tot(base)["TotalHC"] / gdp, 2),
    " VSL-A total", round(100 * tot(base)["VSL_A"] / gdp, 2), "\n")
allh <- run(dscope = "ALL")
cat("Add-on: all-hazard deaths (invasive and chemical included), HC total",
    round(tot(allh)["TotalHC"] / 1e9, 2), " VSL-A total", round(tot(allh)["VSL_A"] / 1e9, 2), "\n")

# ---- 2. Monte Carlo -------------------------------------------------
lnW <- function(lo, hi) list(mu = log(sqrt(lo * hi)), s = (log(hi) - log(lo)) / 3.92)  # no rescaling
Li <- lnW(cty$ill_dia_low, cty$ill_dia_high); Ld <- lnW(cty$dth_dia_low, cty$dth_dia_high)
cat("\n   Simulated mean / WHO mean: diarrhoeal illnesses",
    round(sum(exp(Li$mu + Li$s^2 / 2)) / sum(gsum(ig, "DIA", "mean")), 3),
    " diarrhoeal deaths", round(sum(exp(Ld$mu + Ld$s^2 / 2)) / sum(gsum(dg, "DIA", "mean")), 3), "\n")
sim <- function(rho_c = RHO_C, rho_id = RHO_ID, n = N) {
  out <- matrix(NA, n, 9); colnames(out) <- names(tot(base))
  for (i in 1:n) {
    zc1 <- rnorm(1); zc2 <- rho_id * zc1 + sqrt(1 - rho_id^2) * rnorm(1)
    zi <- sqrt(rho_c) * zc1 + sqrt(1 - rho_c) * rnorm(nc)
    zd <- sqrt(rho_c) * zc2 + sqrt(1 - rho_c) * rnorm(nc)
    pc <- rbeta(1, bc["a"], bc["b"]); ph <- min(rbeta(1, bh["a"], bh["b"]), pc * 0.9)
    zo <- rnorm(1); zi2 <- rnorm(1)                 # one draw per setting; components within a setting move together
    mu <- function(k, z) exp(SDU[k] * z - SDU[k]^2 / 2)
    u <- list(dm_o = U$dm_o * mu("dm_o", zo), dnm_o = U$dnm_o * mu("dnm_o", zo),
              dm_i = U$dm_i * mu("dm_i", zi2), dnm_i = U$dnm_i * mu("dnm_i", zi2))
    do <- P$d_out * exp(rnorm(1, 0, SD_DAYS) - SD_DAYS^2 / 2)
    di <- runif(1, 2.0, 6.2)                       # study range, Zimmermann 2019
    ag <- sample(25:55, 1)                         # adult age at death, assumed 25 to 55
    wm <- exp(rnorm(1, 0, SD_W) - SD_W^2 / 2)
    ill_d <- exp(Li$mu + Li$s * zi) / cty$sh_dia   # diarrhoeal draw expressed on the all-illness scale
    dth_d <- exp(Ld$mu + Ld$s * zd) / cty$dsh_dia
    out[i, ] <- tot(cost(ill_d, dth_d, pc, ph, u, do, di, wm, age_ad = ag))
  }
  out
}
mc <- sim()
q <- function(x) round(quantile(x, c(0.5, 0.025, 0.975)) / 1e9, 2)
cat("\n2. MONTE CARLO, median (2.5th, 97.5th percentile), USD billion, N =", N, "\n")
res <- t(apply(mc, 2, q)); colnames(res) <- c("median", "lo95", "hi95"); print(res)
cat("   Mean of HC total:", round(mean(mc[, "TotalHC"]) / 1e9, 2), " Mean of VSL-A total:", round(mean(mc[, "VSL_A"]) / 1e9, 2), "\n")
for (cs in list(c(0, 0.5), c(1, 0.5), c(0.5, 0), c(0.5, 1))) {
  m2 <- sim(rho_c = cs[1], rho_id = cs[2], n = 3000)
  cat("   rho_c =", cs[1], " rho_id =", cs[2], " HC total:", paste(q(m2[, "TotalHC"]), collapse = " / "),
      "  VSL-A total:", paste(q(m2[, "VSL_A"]), collapse = " / "), "\n")
}

# ---- 3. one-way and scenario analysis ------------------------------
U2_MED   <- list(dm_o = 6.64, dnm_o = 2.98, dm_i = 64.27, dnm_i = 14.21)        # Baral Table 2 medians
U2_MEAN  <- list(dm_o = 12.54, dnm_o = 4.77, dm_i = 117.58, dnm_i = 16.99)      # Baral Table 2 means
# West African study values: outpatient USD 60 per uncomplicated case (Sanni et al. 2024a, 2020 USD, inflated to 2023);
# inpatient USD 113.15, midpoint of 91.13 and 135.16 per admission (Wabada & Abubakar 2020, price year as reported)
U_WA <- list(dm_o = 60 * cpi(2023) / cpi(2020) / INFL_2015_2023, dnm_o = 0,
             dm_i = mean(c(91.13, 135.16)) / INFL_2015_2023, dnm_i = 0)
sc <- function(nm, ...) { t <- tot(run(...)) / 1e9
  data.frame(scenario = nm, HC = t["TotalHC"], VSL_A = t["VSL_A"], VSL_B = t["VSL_B"], VSL_C = t["VSL_C"],
             VSL_D = t["VSL_D"], row.names = NULL) }
scen <- rbind(
  sc("Base case"),
  sc("All-hazard deaths valued (invasive and chemical added)", dscope = "ALL"),
  sc("Invasive illnesses costed as diarrhoea", inf = TRUE),
  sc("Output per worker halved (labour share proxy)", wmult = 0.5),
  sc("Values at 2021 GDP per capita (burden year)", gdppc = cty$gdppc_usd_2021),
  sc("Human capital mortality net of household consumption", net = TRUE),
  sc("Hospitalised share 1.3% (lower CI)", p_hosp = 0.013),
  sc("Hospitalised share 6.9% (upper CI)", p_hosp = 0.069),
  sc("Care seeking 62% (Gambia DHS)", p_care = 0.62),
  sc("Care seeking 33% (Wiens 2025 pooled)", p_care = 0.33),
  sc("Inpatient stay 2.0 days", d_in = 2.0),
  sc("Inpatient stay 6.2 days", d_in = 6.2),
  sc("Outpatient illness 5 days (Saha 2013, The Gambia)", d_out = 5),
  sc("Unit costs at lowest published values (Baral Table 2)", u = list(dm_o = 0.92, dnm_o = 0.47, dm_i = 24.24, dnm_i = 1.18)),
  sc("Unit costs at highest published values (Baral Table 2)", u = list(dm_o = 58.23, dnm_o = 18.87, dm_i = 341.27, dnm_i = 39.75)),
  sc("Unit costs, Baral Table 2 medians", u = U2_MED),
  sc("Unit costs, Baral Table 2 means", u = U2_MEAN),
  sc("Unit costs, West African study values", u = U_WA),
  sc("Adult age at death 30", age_ad = 30),
  sc("Adult age at death 50", age_ad = 50),
  sc("Adult age at death 55", age_ad = 55),
  sc("Discount rate 0%", r = 0),
  sc("Discount rate 5%", r = 0.05),
  sc("Discount rate 7%", r = 0.07),
  sc("Caregiver time for under-5 cases added", caregiver = TRUE),
  sc("All countries at WHO lower bounds simultaneously", ill_all = cty$ill_low, dth_all = cty$dth_low),
  sc("All countries at WHO upper bounds simultaneously", ill_all = cty$ill_high, dth_all = cty$dth_high))
cat("\n3. SCENARIOS (USD billion)\n"); print(scen, digits = 3, row.names = FALSE)

# ---- 4. country results and benchmarks ------------------------------
cc <- sapply(base[c("DMC", "DNMC", "Prod", "MortHC", "MortA")], identity) / 1e6
rownames(cc) <- cty$iso3
cc <- cbind(cc, TotalHC = rowSums(cc[, 1:4]), TotalVSL = rowSums(cc[, c(1:3, 5)]),
            pct_GDP_HC = 100 * rowSums(cc[, 1:4]) * 1e6 / cty$gdp,
            pct_GDP_VSL = 100 * rowSums(cc[, c(1:3, 5)]) * 1e6 / cty$gdp,
            perCap_HC = rowSums(cc[, 1:4]) * 1e6 / cty$pop)
cat("\n4. COUNTRY RESULTS (USD million; % of GDP; USD per person)\n"); print(round(cc[order(-cc[, "TotalHC"]), ], 2))
cat("\n   Benchmarks from included studies:\n")
cat("   Burkina Faso, this model HC total: USD", round(cc["BFA", "TotalHC"]), "m; in int$:",
    round(cc["BFA", "TotalHC"] / cty$pr[cty$iso3 == "BFA"]),
    "m; van Wagenberg & Havelaar (2023), 3 bacterial hazards, 2017: Int$391 m\n")
cat("   Nigeria, this model: productivity USD", round(cc["NGA", "Prod"]), "m; HC mortality", round(cc["NGA", "MortHC"]),
    "m; VSL-A mortality (USD)", round(cc["NGA", "MortA"]), "m; Grace et al. (2018): productivity 5,350 and VSL mortality 16,000 (USD m)\n")
bn <- run(dscope = "ALL"); i_n <- which(cty$iso3 == "NGA")
cat("   Nigeria, all-hazard deaths: HC mortality", round(bn$MortHC[i_n] / 1e6), "m; VSL-A mortality", round(bn$MortA[i_n] / 1e6), "m\n")
cat("   Under-5 share of deaths:", round(sum(cty$dth_u5) / sum(cty$dth_mean), 3), "\n")
cat("   VSL-A per death, int$, range:", paste(round(range(cty$vsl_a)), collapse = " to "), "\n")
cat("   GDP per person aged 15+ proxy W, USD/yr, range:", paste(round(range(cty$gdppc_usd / cty$emp_share)), collapse = " to "), "\n")
cat("   Market USD per int$, range:", paste(round(range(cty$pr), 2), collapse = " to "), "\n")
cat("   Nigeria share of HC total (%):", round(100 * cc["NGA", "TotalHC"] / sum(cc[, "TotalHC"]), 1), "\n")
cat("   Regional HC total at 5% discount (bn):", round(tot(run(r = 0.05))["TotalHC"] / 1e9, 2), "\n")

# ---- 5. hypothesis tests (across the 15 countries, descriptive of the model) ----
cat("\n5. TESTS\n")
grid <- NULL
for (r in c(0, 0.03, 0.05)) for (v in c("A", "B", "C", "D")) {
  x <- run(r = r); mm <- x$DMC + x$DNMC + x$Prod
  sh_hc <- x$MortHC / (mm + x$MortHC); sh_v <- x[[paste0("Mort", v)]] / (mm + x[[paste0("Mort", v)]])
  w <- wilcox.test(sh_hc, sh_v, paired = TRUE, exact = TRUE, conf.int = TRUE)
  grid <- rbind(grid, data.frame(discount = r, vsl = v, med_share_HC = median(sh_hc), med_share_VSL = median(sh_v),
                                 V = unname(w$statistic), p = w$p.value, HL_shift = unname(w$estimate),
                                 n_VSL_higher = sum(sh_v > sh_hc)))
}
cat("   H2 grid, mortality share of total cost, HC vs VSL (paired Wilcoxon, n = 15)\n"); print(grid, digits = 3, row.names = FALSE)
sh_hc <- cc[, "MortHC"] / cc[, "TotalHC"]; sh_vsl <- cc[, "MortA"] / cc[, "TotalVSL"]
cat("   Base H2: HC median", round(median(sh_hc), 3), " VSL-A median", round(median(sh_vsl), 3), "\n")
print(wilcox.test(sh_hc, sh_vsl, paired = TRUE, exact = TRUE, conf.int = TRUE))
cat("   H3: GDP per capita vs cost as % of GDP\n")
h3 <- function(nm, y) { t <- cor.test(cty$gdppc_usd, y, method = "spearman", exact = FALSE)
  cat("   ", nm, ": rho =", round(t$estimate, 3), " p =", signif(t$p.value, 3), "\n"); data.frame(measure = nm, rho = unname(t$estimate), p = t$p.value) }
h3res <- rbind(h3("HC total", cc[, "pct_GDP_HC"]), h3("VSL-A total", cc[, "pct_GDP_VSL"]),
               h3("Productivity + HC mortality only", 100 * (base$Prod + base$MortHC) / cty$gdp),
               h3("HC total, West African unit costs",
                  { x <- run(u = U_WA); 100 * (x$DMC + x$DNMC + x$Prod + x$MortHC) / cty$gdp }))
write.csv(cbind(country = cty$iso3, round(cc, 3)), "tables/country_results.csv", row.names = FALSE)
write.csv(scen, "tables/scenarios.csv", row.names = FALSE)
write.csv(res, "tables/monte_carlo_summary.csv")
write.csv(grid, "tables/h2_grid.csv", row.names = FALSE)
write.csv(h3res, "tables/h3_results.csv", row.names = FALSE)

# ---- 6. reconciliation with the included studies (model values for the same country) ----
ib <- which(cty$iso3 == "BFA"); ing <- which(cty$iso3 == "NGA")
xd <- run(); xa <- run(dscope = "ALL")
m <- function(x, k, i) x[[k]][i] / 1e6
rec <- data.frame(
  item = c("Burkina Faso total, HC (int$ m)", "Burkina Faso total, VSL-A (int$ m)", "Burkina Faso productivity (int$ m)",
           "Nigeria productivity (USD m)", "Nigeria VSL mortality (USD m), diarrhoeal deaths", "Nigeria VSL mortality (USD m), all-hazard deaths",
           "Outpatient direct cost per case (2023 USD), Nigeria", "Inpatient direct cost per admission (2023 USD), Nigeria"),
  model = c((m(xd, "DMC", ib) + m(xd, "DNMC", ib) + m(xd, "Prod", ib) + m(xd, "MortHC", ib)) / cty$pr[ib],
            (m(xd, "DMC", ib) + m(xd, "DNMC", ib) + m(xd, "Prod", ib) + m(xd, "MortA", ib)) / cty$pr[ib],
            m(xd, "Prod", ib) / cty$pr[ib], m(xd, "Prod", ing), m(xd, "MortA", ing), m(xa, "MortA", ing),
            (U$dm_o[ing] + U$dnm_o[ing]) * INFL_2015_2023, (U$dm_i[ing] + U$dnm_i[ing]) * INFL_2015_2023),
  study = c(391, 391, 275, 5350, 16000, 16000, 60 * cpi(2023) / cpi(2020), mean(c(91.13, 135.16))))
rec$ratio <- rec$model / rec$study
cat("\n6. RECONCILIATION WITH INCLUDED STUDIES\n"); print(rec, digits = 3, row.names = FALSE)
cat("   Burkina Faso, share of total from productivity (model HC):",
    round(m(xd, "Prod", ib) / (m(xd, "DMC", ib) + m(xd, "DNMC", ib) + m(xd, "Prod", ib) + m(xd, "MortHC", ib)), 3), "\n")
cat("   Productivity as % of regional GDP:", round(100 * sum(xd$Prod) / gdp, 3), "\n")
write.csv(rec, "tables/reconciliation.csv", row.names = FALSE)
