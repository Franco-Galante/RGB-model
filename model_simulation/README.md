# **RGB Model Simulation**

The C code `gptsim.c` simulates the dynamics of information flow and quality within a simplified ecosystem involving human intervention, training set, the World Wide Web, a Generative AI model, and a search engine, to what concerns a novel topic. It tracks the movement and characteristics of "balls" representing information with quality attributes (RGB components).

## Installation

The simulator requires
- C compiler (`gcc`).
- GNU Scientific Library (GSL) for random number generation and statistical analysis.

## Compilation

To compile, ensure you have the GSL library installed and use:

```bash
gcc gptsim.c -lm -lgsl -lgslcblas -o gptsim -Wall -O3
```

To run the simulator simply:

```bash
./gptsim
```
