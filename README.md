# Master Thesis

*Multiple imputation for the longitudinal networks: a sensitivity study*

## Plan 

1. Find the datasets and estimate the models (same effects)
   - 3 waves, 30 actors, 2 attributes
   - 3 waves, 60 actors, 2 attributes
2. Leave the 1st wave as it is and generate 100 networks for the 2nd and 3rd waves, for each dataset
3. Obtain 2 datasets (100x3 each) and generate different amounts of missing data - low (10%), high (20%), very high (30%)
4. Impute each dataset 4 times (= number of hypothesis)
5. Combine the results.
6. Compare the imputed models to the true es

So it is ((100 nets X 3 missings options X 4 hypothesis) + 2 complete models) X 2 datasets = **2404 models.**

**OR** if we use only 2 missings options, it is **1604 models**.

### Hypothesis

4 models, combining selection no/yes and influence no/yes, where the "yes" is with a parameter value so that power is high but not extreme, e.g., about 80% power.



