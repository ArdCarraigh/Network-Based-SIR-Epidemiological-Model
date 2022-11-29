## WHAT IS IT?

Network-Based SIR Epidemiological Model

## HOW IT WORKS

A first random agent gets infected, then the disease is transmitted depending on social connections linking the agents. Infected agents are then being recovered after a number of days set in the interface. A day is passed when the code has looped through each individual to decide to infect it or not during that day.

## HOW TO USE IT

Parametrize the model interface, click on setup and go.

Parameters are defined as follows:

- simulations-number: Set the number of simulations. Int. Range = [0:10000] by default, but should be adapted to the usage.
- adhesion-limit?: Enable/disable the use of the adhesion limit. Bool.
- time-of-adhesion-limit: Set the adhesion limit. When the number of individuals avoiding infection in a row is superior or equal to the adhesion limit, the code stops. Int. Range = [0:100] by default, but should be adapted to the usage.
- mimetism?: Enable/Disable the use of the number of infected individuals in computing probablity of transmission. Bool.
- aff?: Enable/Disable the use of the social network links in computing probability of transmission. Bool.
- R0: Set the basic reproductive rate of the pathogen. Float.
- recovery-time: Set the infectious period of the pathogen. Int.
- attributes_: Name of the attributes table to input to the model. String.
- links_: Name of the links table to input to the model. String.

The attributes table consists of 2 columns WITHOUT header: the individual ID and the susceptibility of individuals. The susceptibility is integrated as the size of the nodes, and should be equal to 1 by default. Individual susceptibility might vary between individuals and should be equal to 0 if the insividual is vaccinated.

The links table is an edgelist of social connections between individuals. It consists of 3 columns WITHOUT header: the ID of the individual emitting the link, the ID of the individual receiving it and the strength of the link.

Examples for both tables are provided with the model files.

On the right side of the interface, you have an output with for each simulation in chronological order:

- id of the infected individual
- time spent since last infection
- time spent since the beginning of the simulation
- the chronological rank of infection
- the day of infection

This output is written in your working directory as "Output.csv" when the code stops running (i.e. when all simulations requested have finished running).

Functions for R Programming Language are provided with the model files in order to ease import and analysis.

On the far right side of the interface you have plots to track the number of individuals in each compartments (Susceptible - Infected - Recovered) as well as the probability of infection of selected individual. The user might edit the plot of probability of infection to track the desired individuals.


## THINGS TO NOTICE

Disabling the visual update allows the model to run faster.
The speed of transmission is relative to the R0 max set in the model. By default it is set to 10 and assumes a probability of transmission of 1 when R0 = 1. The default value moght be adapted to the usage.

## THINGS TO TRY

The user is encouraged to fiddle with R0 and recovery time values in order to observe their respective effects.

## EXTENDING THE MODEL

- Create an infection process that doesn't rely on a arbitrary R0 max.
- Implement a vaccination option, with various strategies (random, trait-based, centrality-based)

## NETLOGO FEATURES

We used netlogo built-in linking functionality to represent the social connetions between individuals.

## RELATED MODELS

Other epidemiological models in the netlogo's model library:

- Biology/Disease Solo
- Biology/HIV
- Biology/Virus 

## CREDITS AND REFERENCES

Model developed by Cédric Sueur, Valéria Romano and Maxime Pierron.

[![DOI](https://zenodo.org/badge/545447277.svg)](https://zenodo.org/badge/latestdoi/545447277)

URL: https://github.com/ArdCarraigh/Network-Based-SIR-Epidemiological-Model
