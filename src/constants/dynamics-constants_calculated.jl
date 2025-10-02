export BETA_MEAN,
    EPSILON

const BETA_MEAN = calculate_beta(R0, SIGMA, GAMMA, MU, POPULATION_N)
const EPSILON = calculate_import_rate(MU, R0, POPULATION_N)
