# Master Thesis

*Multiple imputation for the longitudinal networks: a sensitivity study*

## Plan 

1. Find the datasets and estimate the models (same effects)
   - 3 waves, 30 actors, 2 attributes
   - 3 waves, 60 actors, 2 attributes
2. Leave the 1st wave as it is and generate 100 networks for the 2nd and 3rd waves, for each dataset
3. Obtain 2 datasets and generate different amounts of missing data - low (10%), high (20%), very high (30%) x 3 hypothesis (see below)
4. Impute each dataset 4 times (see thetas)
5. Combine the results
6. Compare the imputed models to the true estimates

So it is (100 nets X 3 missings options X 3 hypothesis) X 2 datasets = **1800 missing data models.**

**OR** if we use only 2 missings options, it is **1200 missing data models.**

Each of then should be imputed 4 times with different thetas:
- 1800 x 4 = 7200 or
- 1200 x 4 = 4800

Possible datasets:
'Teenage Friends and Lifestyle Study' data - select 30 actors and **another** 60 actors from 129 non-missing actors.

### Hypothesis

1. Missings depend on the network
2. Missings depend on the behavior
3. Missings depend on both 

### Thetas

4 models, combining selection no/yes and influence no/yes, where the "yes" is with a parameter value so that power is high but not extreme, e.g., about 80% power.

All selection effects ('altX', 'altSqX', 'egoX', 'egoSqX', 'egoXaltX', 'simX', 'diffX', 'diffSqX', 'higher', 'sameX', 'egoDiffX', 'egoPlusAltX')

All influence effects ('linear', 'quad', 'avAlt', 'avSim', 'totAlt', 'totSim')



