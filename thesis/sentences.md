# Master Shittens
***
## Optimization (NSGA-II)
- If we can find a set of solutions that they don’t dominate each other and not dominated by any other solutions, we call them “Pareto-optimal” solutions.
- Non-dominated Sorting classifies all individuals to different Pareto-optimal front level
- Crowding Distance (explicit diversity preserving mechanism)
-- perform Crowding-sort that uses crowding distance that is related with the density of solutions around each solution. The less dense are preferred
- Domination
A solution x(1) is said to dominate the other solution x(2) if both condition 1 and 2 below are true:
-- Condition 1: x(1) is no worse than x(2) for all objectives
-- Condition 2: x(1) is strictly better than x(2) in at least one objective
- Elitist principle

[ref] Optimization for Engineering Design Algorithms and Examples
- With the advent of computers, optimization has become a part of computer-aided design activities
- not only to achieve just a feasible design, but also a design objective
- optimization algorithms begin with one or more design solutions supplied by the user and then iteratively check new design solutions in order to achieve the true optimum solution
- The purpose of the formulation procedure is to create a mathematical model of the optimal design problem, which then can be solved using an optimization algorithm
- optimization algorithms are routinely used in aerospace design activities to minimize the overall weight, simply because every element or component adds to the overall weight of the aircraft
- a naive optimal design is achieved by comparing a few (limited up to ten or so) alternative design solutions created by using a priori problem knowledge
- efficiency and speed of optimization algorithms depend, to a large extent, on the number of chosen design variables

***
## Turbine Design
### Partial admission

[ref](https://www.conceptsnrec.com/blog/preliminary-sizing-of-supersonic-turbines-with-partial-admission-for-best-performance-using-axial)
- If the blade heights of the various impulse turbines are still too low to deliver acceptable efficiency levels, or be viable for manufacturing, a partial admission is the next step in the consideration of the turbine design
- “Exchange” a reduction in the admission arc for an increase in blade height
- Drawbacks: various additional losses: fill-in, empty-out, unsteadiness in the active jet domains, ventilation/recirculation (or “pumping”), and leakages in active and passive arc domains
- Fully expanded supersonic flow for each bladed passage
- The design target for a supersonic nozzle is to find **throat/exit area ratio** of the passage that delivers minimum loss for the intended velocity triangles of the designed stage
- constraints: aerodynamic design (allowed loadings, diffusion, reactions), structural design (blade and disk stresses) and manufacturability (thicknesses, blade separation distances, blade heights, blade angles)
- Stage input parameter for optmization: mass flow, inlet total temperature and total pressure, and pressure ratio
- Minimum allowed blade height and separation, tip clearance in rotor, no blade *twist* (constant blade profile hub to tip)
- Output: admission level,throat to exit area ratio for a given Mach at exit
- Radial sizing and rotation speed of the stage are adjusted for best performance
- Mach number at exit must coincide with Mach number for fully expanded supersonic flow, which is an indication of correct setuo of throat/exit area ratio
- Exit swirl angle might be negtive (zero swirl is desirable for better total-to-static efficiency) --> compromise between increasing exhaust energy loss and reduction in vetilation loss
- After optimization, the supersonic stage can be studied for **off-design** performance
- Performance map are plots of efficiencies for several PR at a given speed
