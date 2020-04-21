# Getting Started

## Install dependencies
Use Anaconda to set up an environment: `conda env create -f environment.yml`

Activate the environment: `conda activate covid-group-tests`

## Estimating prevalence
Prevalence estimation takes as input the viral load (copies per ml) measured in a small number of batches (or pools). Together with these, we use Kernel Density Estimates (KDEs) of the viral load distribution observed in a population to build a Maximum Likelihood Estimate of prevalence in the sampled population. In the examples below, we'll use pre-trained KDEs found in the [examples](../../examples/) folder. (These were trained using `train_prevalence_kde.py` on our simulated population.)

### Prevalence example 1
In this example we'll suppose we collected **N=9,216** samples and split them into **b=24** pools each with **n=384** samples. Let's suppose 1/9,216 samples was positive (**p=0.0109%**), with a viral load of ~10,000. In that case, our pooled measurements might look like [example_1.txt](../../examples/observations.viral_load.example_1.txt) (*i.e.*, 23 zeros and 1 measurement with a viral load ~10,000/384).

We run prevalence estimation as follows:

`python estimate_prevalence.py --batch-size 384 --num-batches 24 --observations ../../examples/observations.viral_load.example_1.txt --kde-path ../../examples/KDEs/seir_kde`

The expected output is:

`batch_size      384

num_batches     24

observations    ../../examples/observations.viral_load.example_1.txt

kde_path        ../../examples/KDEs/seir_kde

p0      0.005

Estimated prevalence: 0.0111%`

(Note: p0 is the initial estimate of the EM algorithm. In this case, 0.5%)

### Prevalence example 2
Now let's suppose we collected **N=2,304** samples and split them into **b=24** pools each with **n=96** samples. If 3/2,304 samples were positive (**p=0.1302%**), then our pooled measurements might look like [example_2.txt](../../examples/observations.viral_load.example_2.txt).

Run prevalence estimation:

`python estimate_prevalence.py --batch-size 96 --num-batches 24 --observations ../../examples/observations.viral_load.example_2.txt --kde-path ../../examples/KDEs/seir_kde`

The expected output is:

`batch_size      96

num_batches     24

observations    ../../examples/observations.viral_load.example_2.txt

kde_path        ../../examples/KDEs/seir_kde

p0      0.005

Estimated prevalence: 0.1406%`


