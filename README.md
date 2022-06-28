# Project overview

## What are you trying to do? Articulate your objectives using words that would be familiar to someone who had taken an undergraduate class in the topic.

Adaptive immune receptors are created via a stochastic process called V(D)J recombination. We are trying to verify whether microhomology guides V(D)J recombination events such as trimming and N nucleotide insertions.

## What are you curious about? What will keep pushing you forward when things get tough?


We're curious about how V(D)J recombination works!
Specifically, we're interested if microhomology is important for V(D)J recombination and if so, how many nucleotides of microhomology are optimal for the process.
We are also interested in the relationship between microhomology and trimming/N-insertions as well as the order in which they happen.
Since current V(D)J recombination models are not aware of microhomology, we want to incorporate microhomology information to improve model accuracy.

## How is it done today, and what are the limits of current practice?

1. The current V(D)J recombination models are not aware of microhomology. We expect that incorporating microhomology into model formulation will have a better fit to the data.
2. Currently, trimming and N insertion are assumed to be independent. If microhomology has an influence on the V(D)J recombination events, this will change the way these events are thought of.
3. Although V,D,J gene usage bias has been observed and reported (and learned by the models from point 1), current models learn a V/D/J gene usage bias. We suspect that V/D/J gene usage is at least in part influenced by microhomology between the two joining segments. Therefore, a V/D/J bias might not need to be learned.


## What's new in your approach and why do you think it will be successful?

We will:
1. Formulate the microhomology-dependent V(D)J recombination model as a probabilistic model
2. We will formulate our model based on our understanding of real biological mechanisms


## Who cares? If you're successful, what difference will it make?

V(D)J recombination is at the heart of immune repertoire creation.
Having a better model of how this process happens is important for:

* assessing diversity in individual immune repertoires
* comparing between individual repertoires
* people interested in antibody and TCR sequence generation
* a better understanding of the implications of microhomology in trimming and N insertions. Currently, immunologists think N insertions are mecanism in place to "add diversity" to the repertoire. If it's a mecanism to ensure better end joining via microhomology, this changes the paradigm.


## What are the risks and the payoffs?

There is a risk that we cannot obtain direct evidence for what we look for in the data. However, securing collaborations for in vitro validation could be an option.
There is also a risk of a mixed signal, where microhomology contributes but is not the primary cause. This risk is mitigated by preliminary results we obtained.

The payoffs are in the applications above.


## Are there any other categorically different approaches that could be applied here?

I think our approach is the most basic one.
We could perhaps only have a data-driven approach, with no modelling but it is not categorically different from ours.
??


## Is it possible that better data would make this project irrelevant?

* Immune receptor repertoire data are snapshots of the result after V(D)J recombination. No type of data can currently show the underlying intermediary steps.


## Is there pre-existing work/code that could be leveraged to at least get a first-pass answer for the underlying scientific question? To validate the approach?

* Maggie's work for her trimming model: ??
* Maggie's work in the GWAS paper: ??

Maggie's trimming model contains some important features we can look for wrt microhomology.
The GWAS data is a cool application to look for differences between individuals once we have a model.


## What data will you use? Are there appropriate hold-out sets?

We will use simulated IGoR sequences.
We will use preterm and cord blood repertoires to limit N insertions
We will use thymic TCR data from young children, for a bigger amount of non-productive sequences.
We will use some Omenn's Syndrome repertoires, where they have oligoclonal TCR repertoires (defective VDJ recombination, will keep some of the model's "moving parts" stable)
We will test on some BCR data to verify if the microhomology model applies accross all immune receptors.
We will use the GWAS-Emerson set, to validate in a large cohort and look for genomic variation associations


## Sketch the approach, broken down into steps, with expected amounts of time and intermediate steps for each.

### Preliminary signal check (1 months):

Make sure there is signal in the data. 
* See if we see a signal for less trimming/N insertions in junctions with microhomology
* See if there is a bias in gene usage vs microhomology
* Create the trimming tables:

  trim V
  01234567
t0X
r1
i2
m3
 4
D5
 6
 7
Where X is the number of microhomologous nucleotides once the sequence is trimmed as per i,j coordinates in the table.


#### Test

Do we see signal?

#### How do we decide to move onto the next stage?

If there is signal in the data, sounds reasonable to continue.


### Level 1 model (2 months):

Model only trimming, no N insertions.


#### Test

Does this model fit better the no N insertions data than the current existing model?

#### How do we decide to move onto the next stage?

Depending on the test result, we can move onto modeling sequences with N insertions.


### Level 2 model (2 months):

Add the N insertions component.


#### Test

Does this model fit better the no N insertions data than the current existing model and the model Level 1?

#### How do we decide to move onto the next stage?

Depending on the test result, we can move onto modeling sequences with N insertions.


### Poke at the model (2 months):

Look for someone who is willing to make this model fast.
Poke at the model.
Test on different datasets.
Write manuscript draft.

#### Test
At this stage, if some things need to be verified in vitro, make a concrete plan of what needs to be verified and contact relevant teams for collaboration.
