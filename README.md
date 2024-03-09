## Chromosome Simulation

This project is a genetic algorithms simulation implemented in Python. Genetic algorithms are heuristic search algorithms inspired by the process of natural selection and genetics.

### Overview

The project simulates a genetic algorithm to optimize a mathematical function represented by coefficients \(a\), \(b\), and \(c\). The genetic algorithm aims to find the maximum value of this function within a specified range by evolving a population of binary chromosomes.

### Features

- **Chromosome Initialization** - random binary strings representing chromosomes are generated to form an initial population.
- **Fitness Evaluation** - the fitness of each chromosome is calculated based on the mathematical function with given coefficients.
- **Selection** - chromosomes are selected for reproduction using a fitness-proportionate selection mechanism.
- **Crossover** - selected chromosomes undergo crossover to produce new offspring.
- **Mutation** - a mutation operation introduces random changes in the chromosomes to maintain genetic diversity.
- **Elitism** - the fittest chromosome from the previous generation is preserved in the current population.
- **Evolution** - the algorithm iterates through multiple generations, optimizing the population towards the maximum fitness value.

### How to Use

1. **Input Configuration**: Provide input parameters in the `input.txt` file. Parameters include population size, coefficient values, precision, crossover probability, mutation probability, and number of iterations.
2. **Run the Code**: Execute the Python script (`main.py`), which reads the input parameters, performs the genetic algorithm simulation, and writes the evolution progress to the `evolution.txt` file.
3. **Output Analysis**: Analyze the `evolution.txt` file to observe the evolution of the population over iterations and determine the maximum fitness achieved.
