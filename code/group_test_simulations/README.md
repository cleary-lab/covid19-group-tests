# Getting Started

## Install dependencies
Use Anaconda to set up an environment: `conda env create -f environment.yml`

Activate the environment: `conda activate covid-group-tests`

## Estimating prevalence
Prevalence estimation takes as input the Ct values measured in a small number of batches (or pools). Together with these, we use Kernel Density Estimates (KDEs) of Ct values that follow from the viral load distribution observed in a population to build a Maximum Likelihood Estimate of prevalence in the sampled population. In the examples below, we'll use pre-trained KDEs found in the [examples](../../examples/) folder. (These were trained using `train_prevalence_kde.py` on our simulated population.)

### Prevalence example 1
In this example we'll suppose we collected **N=9,216** samples and split them into **b=24** pools each with **n=384** samples. Let's suppose 1/9,216 samples was positive (**p=0.0109%**), with a viral load of ~10,000. In that case, the Ct values from pooled measurements might look like [example_1.txt](../../examples/observations.Ct_values.example_1.txt) (*i.e.*, 23 pools with undetected RNA and 1 with a Ct value corresponding to a viral load ~10,000/384).

We run prevalence estimation as follows:

```bash
python estimate_prevalence.py --batch-size 384 --num-batches 24 --observations ../../examples/observations.Ct_values.example_1.txt --kde-path ../../examples/KDEs/seir_kde
```

The expected output is:

```bash
batch_size      384
num_batches     24
observations    ../../examples/observations.Ct_values.example_1.txt
kde_path        ../../examples/KDEs/seir_kde
p0      0.005
fp_rate 0.002
CI      False
Estimated prevalence: 0.0105%
```

(Note: p0 is the initial estimate of the EM algorithm. In this case, 0.5%)

You can also use the --CI flag to compute bootstrap confidence intervals, though this substantially increases runtime.

### Prevalence example 2
Now let's suppose we collected **N=2,304** samples and split them into **b=24** pools each with **n=96** samples. If 3/2,304 samples were positive (**p=0.1302%**), then our pooled measurements might look like [example_2.txt](../../examples/observations.Ct_values.example_2.txt).

Run prevalence estimation:

```bash
python estimate_prevalence.py --batch-size 96 --num-batches 24 --observations ../../examples/observations.Ct_values.example_2.txt --kde-path ../../examples/KDEs/seir_kde
```

The expected output is:

```bash
batch_size      96
num_batches     24
observations    ../../examples/observations.Ct_values.example_2.txt
kde_path        ../../examples/KDEs/seir_kde
p0      0.005
fp_rate 0.002
CI      False
Estimated prevalence: 0.1360%
```

## Individual case identification
Our case identification algorithm takes as input the positive or negative result of each pooled test, along with the composition of each pool. The output is a list of putatively positive samples to be validated in a second stage of testing.

### Identification Example 1
Here we'll use the [pools](../../examples/pool_composition.example_3.txt) and [results](../../examples/observations.positive_negative.example_3.txt) of real tests we performed on de-identified samples. This experiment consisted of 48 samples -- the first 47 samples were negative and the **48th sample was positive**. Each sample was split into 3 out of 6 total pools.

Run the code for sample identification as follows:

```bash
python identify_positive_samples.py --observations ../../examples/observations.positive_negative.example_3.txt --pools ../../examples/pool_composition.example_3.txt
```

The expected output is:

```bash
observations    ../../examples/observations.positive_negative.example_3.txt
pools   ../../examples/pool_composition.example_3.txt

Found 2 putative positives:
sample 32
sample 48
```

Each of these samples was included in all 3 positive pools (1, 2, and 4).

