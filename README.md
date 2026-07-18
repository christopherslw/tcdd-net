# Estimation of mobility trajectory-constrained dynamic dependence networks

Implements the following procedure to estimate spatial dependence between locations or points of interest
from raw mobility trajectories:

Step 1: Compute a pilot visitation graph by counting visits from each location.

Step 2: Fit an autoregressive model and obtain residuals for each location.

Step 3: Select edges from the pilot fit via group lasso regression on the residuals, grouped across all time lags.
