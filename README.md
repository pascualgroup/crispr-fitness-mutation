# Stochastic CRISPR model implementation

## Running a single simulation

There are three categories of parameters, organized into separate structs in the code:

1. Run parameters (`RunParameters` struct): start and stop time of the simulation, and random seed (`rng_seed`)
2. Initialization parameters (`InitializationParameters` struct): parameters governing the initial population
3. Model parameters (`Parameters` struct): parameters governing the dynamical process

You can run a single simulation using `run.jl`, with a single argument giving the location of the parameters file in JSON format, like this:

```sh
julia /path/to/stochastic-crispr/julia/main.jl parameters.json
```

An example parameters file is in `example/parameters.json`.

## Running a parameter sweep

To sweep across many parameter values, you'll want to write a script to generate a working directory for each simulation, containing a parameters file for that simulation and perhaps also a script that submits the simulation to a cluster.

An example of this script is `/simulation/generate-sweep.jl`, which calls a JSON file designating the parameter intervals (`sweep-parameters.json`). This script is written for a cluster operating with SLURM. Number of processes, cores and jobs can be customized. Run
```sh
julia /simulation/generate-sweep.jl
```

## Model structure and simulation algorithm

The model is a stochastic, continuous-time, discrete-event simulation (Gillespie-style) that steps from event to event, with a random amount of time, real-valued, between events.

## Model state

Model state is represented by the `State` struct, which contains two fields, `bstrains`, of type `BStrains`, and `vstrains`, of type `VStrains`.

The `BStrains` and `VStrains` structs are very similar.
Conceptually, they represent a collection of microbial/viral strains, and although it would make conceptual sense to just have `bstrains::Vector{BStrain}`, for efficiency of memory layout they are implemented using parallel arrays.

## Initialization

## Event loop

Each step follows this algorithm:

1. The time to the next event is computed by drawing a waiting time from a Poisson process whose rate is the total rate of all possible events in the system.
2. The type of event is chosen proportional to the total rate of all events of that type (e.g., all possible microbial growth events).
3. For each event type, the specific event is chosen by sampling proportional to the rates of the individual events. For example, for microbial growth events, since microbial cells grow at the same rate, the growth rate of particular strain is proportional to the abundance of that strain, and therefore a strain is chosen proportional to its abundance.
4. The event is executed.

There are four events:

1. Microbial growth: `const MICROBIAL_GROWTH = 1`
2. Microbial washout ("death"): `const MICROBIAL_DEATH = 2`
3. Viral decay: `const VIRAL_DECAY = 3`
4. Contact: `const CONTACT = 4`
4. MICROBIAL immigration: `const MICROBIAL_IMMIGRATION = 5`

### Microbial growth

1. Randomly select a microbial strain with probability proportional to growth rate

### Contact Event

1. Randomly select a bacterial strain and a viral strain proportional to population size
2. If immune, infect with probability p_crispr_failure_prob
3. If not immune, defend successfully & acquire spacer with prob q_spacer_acquisition_prob

Infect:

1. Draw the number of mutations for each of `beta` newly created virus particles; the number is `Binomial(n_pspacers, mu)`.
2. Unmutated virus particles just increase the abundance of the existing viral strain.
3. Mutated virus particles acquire newly created protospacers at randomly selected loci.

## Generating Figures

Please request simulation and analysis ouputs from author, or generate such output on a cluster. Note that such files range from 50 to 200 GB each. Therefore storage resources are crucial to be allocated beforehand. For assistance email armun@uchicago.edu. 

Outputs will be organized into databases (.sqlite) that should be placed in the same directory as the following scripts listed below. 

### Figure 1, 2 & 3 and Supplementary Figure 1 & 2
Run the following
```sh
python figures_1-2-3_and_supp-figures_1-2.py
```

### Figure 4 and Supplementary Figures 3 & 4
Run the following
```sh
python figure_4_and_supp-figures_3-4.py
```

### Figure 5 and Supplementary Figures 5-7
Run the following
```sh
python figure_5_and_supp-figures_5-7.py
```




