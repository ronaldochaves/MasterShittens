## Master Shittens
# Optimization (NSGA-II)
- If we can find a set of solutions that they don’t dominate each other and not dominated by any other solutions, we call them “Pareto-optimal” solutions.
- Non-dominated Sorting classifies all individuals to different Pareto-optimal front level
- Crowding Distance (explicit diversity preserving mechanism)
-- perform Crowding-sort that uses crowding distance that is related with the density of solutions around each solution. The less dense are preferred
- Domination
A solution x(1) is said to dominate the other solution x(2) if both condition 1 and 2 below are true:
-- Condition 1: x(1) is no worse than x(2) for all objectives
-- Condition 2: x(1) is strictly better than x(2) in at least one objective
- Elitist principle