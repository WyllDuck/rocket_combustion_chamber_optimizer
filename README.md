# rocket_combustion_chamber_optimizer

Calculation of the dimensions of the regenerative cooled combustion chanmber of a rocket engine
Assumptions:
- Methane coolant at supercritical pressure (change of phase does not requiere an extra heat transfer)
- Cooled from the throat of the engine to the injection plate



# Explanaition
This code discritizes the combustion chamber in N section (~20) and solves the temperature at the combustion chamber wall (T_wi), the coolant channel wall (T_wo) and the coolants temepratures (T_in, T_out). The concept to solve this system was proven in an excel sheet for 1 section and now I want to scale it up to get the final design parameters.

The final objective is to determine the temperature along the CC walls to find a feasible design for the CC. 

Questions:
- I have the combustion temperature of the hot gas (inside the combustion chamber) and the inner diameter of the CC (Di). With this information I want to obtain the temperature at every position of the CC (inside the thrust chamber and the nozzle) based on the Di using something similar to the ideal gas law. Could that work? Which function is that one?

