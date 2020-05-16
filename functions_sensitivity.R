require(ggplot2)
require(data.table)
require(forcats)
require(tidyverse)
require(patchwork)
require(deSolve)
require(ggrepel)

## Function to simulate observed prevalence, adjusted for Se's of clinical symptoms (asx/mild/severe)
sim.prev_symptoms <- function(Se.asx, Se.mild, Se.severe, Sp=1, prop.asx, prop.mild, prop.severe, true.prev) {
	
	return(sum((
		(Se.asx*prop.asx)+
		(Se.mild*prop.mild)+
		(Se.severe*prop.severe))*true.prev)+
		((1-Sp)*(1-true.prev)))
	
}

## Function to simulate observed prevalence, adjusted for Se's of times since infection (1/2/3)
sim.prev_TSI <- function(Se.1, Se.2, Se.3, Sp=1, prop.1, prop.2, prop.3, true.prev) {
	
	return(sum((
		(Se.1*prop.1)+
		(Se.2*prop.2)+
		(Se.3*prop.3))*true.prev)+
		((1-Sp)*(1-true.prev)))
	
}

## Function to simulate observed prevalence, adjusted for Se's of clinical symptoms (asx/mild/severe) & times since infection (1/2/3)
sim.prev_TSI_symptoms <- function(Se.asx.1, Se.asx.2, Se.asx.3, Se.mild.1, Se.mild.2, Se.mild.3, Se.severe.1, Se.severe.2, Se.severe.3, Sp=1, prop.asx.1, prop.asx.2, prop.asx.3, prop.mild.1, prop.mild.2, prop.mild.3, prop.severe.1, prop.severe.2, prop.severe.3, true.prev) {
	
	return(sum((
		(Se.asx.1*prop.asx.1)+
		(Se.asx.2*prop.asx.2)+
		(Se.asx.3*prop.asx.3)+
		(Se.mild.1*prop.mild.1)+
		(Se.mild.2*prop.mild.2)+
		(Se.mild.3*prop.mild.3)+
		(Se.severe.1*prop.severe.1)+
		(Se.severe.2*prop.severe.2)+
		(Se.severe.3*prop.severe.3))*true.prev)+
		((1-Sp)*(1-true.prev)))
	
}

## Function to calculate the weighted Se in the population, given the Se's and prevalence of clinical symptoms (asx/mild/severe) in the positive controls
Se.weighted_symptoms <- function(Se.asx, Se.mild, Se.severe, v.asx, v.mild, v.severe) {
	
	return(
		(Se.asx*v.asx)+
		(Se.mild*v.mild)+
		(Se.severe*v.severe))
	
}

## Function to calculate the weighted Se in the population, given the Se's and prevalence of times since infection (1/2/3) in the positive controls
Se.weighted_TSI <- function(Se.1, Se.2, Se.3, v.1, v.2, v.3) {
	
	return(
		(Se.1*v.1)+
		(Se.2*v.2)+
		(Se.3*v.3))
	
}

## Function to calculate the weighted Se in the population, given the Se's and prevalence of clinical symptoms (asx/mild/severe) & times since infection (1/2/3) in the positive controls
Se.weighted_symptoms_TSI <- function(Se.asx.1, Se.asx.2, Se.asx.3, Se.mild.1, Se.mild.2, Se.mild.3, Se.severe.1, Se.severe.2, Se.severe.3, v.asx.1, v.asx.2, v.asx.3, v.mild.1, v.mild.2, v.mild.3, v.severe.1, v.severe.2, v.severe.3) {
	
	return(
		(Se.asx.1*v.asx.1)+
		(Se.asx.2*v.asx.2)+
		(Se.asx.3*v.asx.3)+
		(Se.mild.1*v.mild.1)+
		(Se.mild.2*v.mild.2)+
		(Se.mild.3*v.mild.3)+
		(Se.severe.1*v.severe.1)+
		(Se.severe.2*v.severe.2)+
		(Se.severe.3*v.severe.3))
	
}

## Function to adjust the observed prevalence for a single Se & Sp
adj.func <- function(prev.obs, Se, Sp) {
	
	return((prev.obs+Sp-1)/(Se+Sp-1))
	
}

## Function to model a time step for the SEIR model
dx.dt.SEIR <- function(t, y, params) {
	
	## Calculate the change in Susceptible
	dS <- -params["beta"]*y["S"]*y["I"]
	
	## Calculate the change in Exposed
	dE <- params["beta"]*y["S"]*y["I"] - params["sigma"]*y["E"]

	## Calculate the change in Infectious
	dI <- params["sigma"]*y["E"] - params["gamma"]*y["I"]
	
	## Calculate the change in Recovered
	dR <- params["gamma"]*y["I"]
	
	## Track cumulative incidence
	dC <- params["beta"]*y["S"]*y["I"]
	
	## Get changes in S, E, I, R, and cumulative incidence at the current time step
	return(list(c(dS, dE, dI, dR, dC)))
	
}
