const POPULATION_N = 500_000
const LATENT_PER_DAYS = 10
const DUR_INF_DAYS = 8
const R0 = 16.0
const SIGMA = 1 / LATENT_PER_DAYS
const GAMMA = 1 / DUR_INF_DAYS
const ANNUAL_BIRTHS_PER_K = 27
const LIFE_EXPECTANCY_YEARS = 1000 / ANNUAL_BIRTHS_PER_K
const MU = 1 / (LIFE_EXPECTANCY_YEARS * 365)
const BETA_MEAN = calculate_beta(R0, GAMMA, MU, 1, POPULATION_N)
const BETA_FORCE = 0.0
const EPSILON = calculate_import_rate(MU, R0, POPULATION_N)
const VACCINATION_COVERAGE = 0.8
