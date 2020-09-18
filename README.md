# coverage-path-generation
Designed for Catabot at robotics lab of Dartmouth College. 
This algarithm intends to generate the coverage path for an ASV(Catabot) to measure and take photos of a known, irregular area. Compared with previous methods, this algarithm features the adaptiveness to bathymatry. In other words, the interval distance of the path varies with the depth of the water. The algarithm is achieved by conducting area decomposition and integar programing. Gurobi solver is implemented to solve the integar programing and give the optimal sequence of visiting each decomposed area.
Future work:   
generalize the scenario to multi-agent.
