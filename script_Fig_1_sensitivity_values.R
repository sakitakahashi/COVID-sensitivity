## Read in source code
source("functions_sensitivity.R")

### Set the sensitivity in each category
Se.asx <- 0.4
Se.mild <- 0.6
Se.severe <- 0.95

################################################################################

## Generate scenarios
data.frame(which_scenario="1", prop.asx=0.43, prop.mild=0.52, prop.severe=0.05) %>%
	add_row(which_scenario="2", prop.asx=0.05, prop.mild=0.45, prop.severe=0.50) %>%
	add_row(which_scenario="3", prop.asx=0, prop.mild=0, prop.severe=1) -> scenarios

scenarios$Se_in_controls <- NULL

## Calculate the Se in the positive control sets
for(i in 1:nrow(scenarios)) {
	
	Se.weighted_symptoms(
		Se.asx=Se.asx,
		Se.mild=Se.mild,
		Se.severe=Se.severe,
		v.asx=scenarios$prop.asx[i],
		v.mild=scenarios$prop.mild[i],
		v.severe=scenarios$prop.severe[i]) -> scenarios$Se_in_controls[i]
	
}

scenarios

## Calculate the Se in the general population
Se.weighted_symptoms(
	Se.asx=Se.asx,
	Se.mild=Se.mild,
	Se.severe=Se.severe,
	v.asx=scenarios$prop.asx[1],
	v.mild=scenarios$prop.mild[1],
	v.severe=scenarios$prop.severe[1])
