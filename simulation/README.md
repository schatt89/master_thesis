# Simulations

This readme contains the description of the models used for simulations

All selection effects ('altX', 'altSqX', 'egoX', 'egoSqX', 'egoXaltX', 'simX', 'diffX', 'diffSqX', 'higher', 'sameX', 'egoDiffX', 'egoPlusAltX')

All influence effects ('linear', 'quad', 'avAlt', 'avSim', 'totAlt', 'totSim')

## Simulation1.R

Data: Glasgow dataset (129 non-missing actors)
* Smaller data: first 30 actors
* Bigger data: last 60 actors
* 3 waves
* Covars: gender, tobacco[,1]
* Behavior: alcohol

Model specifications:
* The effects are taken from the SimulateNetworkaBehavior functions
* Structural: transTrip, cycle3
* Covars: egoX, altX, simX
* Behavior: avAlt interaction=("friendship"), avAlt and quad interaction=("", "friendship")

Both models converge.

## Simulation2.R (gonna try)

Data: Glasgow dataset (129 non-missing actors)
* Smaller data: first 30 actors
* Bigger data: last 60 actors
* 3 waves
* Covar: gender
* Behavior: alcohol, tobacco

Model specifications:

* Structural: transTrip, cycle3
* Covars: sameX
* Behavior: [altX, egoX, egoXaltX] for selection, [quad, avAlt, avSim] for influence


