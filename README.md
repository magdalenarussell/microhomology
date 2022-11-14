# Project overview

## What are you trying to do? Articulate your objectives using words that would be familiar to someone who had taken an undergraduate class in the topic.

Adaptive immune receptors are created via a stochastic process called V(D)J recombination. 
During V(D)J recombination, sequences are joined and ligated together through the non-homologous end joining (NHEJ) process.
It is fairly well established that microhomology (MH) between the two sequences being joined occurs for most NHEJ events, except for when direct ligation occurs.
Mechanistically, we are interested in the relationship between V(D)J recombination events (such as trimming and N-insertion) and MH; for example, to what extent does internal MH influence trimming?

## What are you curious about? What will keep pushing you forward when things get tough?

We're curious about how V(D)J recombination works!
Specifically, we're interested in whether microhomology is important for V(D)J recombination and if so, how many nucleotides of microhomology are optimal for the process.
We are also interested in the relationship between microhomology and trimming/N-insertion, as well as the order in which they happen.
Because this relationship has not been fully defined in a human system, any result (either positive or negative) will be informative. 

## How is it done today, and what are the limits of current practice?

Established probabilistic models of V(D)J recombination, and individual V(D)J recombination events (such as nucleotide trimming), are not aware of microhomology.
Currently, trimming and N insertion are assumed to be independent.
If microhomology has an influence on these V(D)J recombination events, it will change the way these events are thought of.
As such, we expect that incorporating microhomology into model formulation may lead to a better fit to the data and inform our understanding of the biological mechanism.

## What's new in your approach and why do you think it will be successful?

We will:
1. Formulate microhomology-dependent models of various V(D)J recombination events as probabilistic models; We will start by modeling trimming (see detailed plan below)
2. We will formulate our model based on our understanding of real biological mechanisms

## Who cares? If you're successful, what difference will it make?

V(D)J recombination is at the heart of immune repertoire creation.
Many people care about antibody and TCR sequence generation.
Having a better model of how this process happens is important for:

* understanding how the V(D)J recombination mechanism works in healthy humans
* understanding the implications of microhomology for trimming and N-insertions. Currently, immunologists think N-insertions are a mechanism in place to "add diversity" to the repertoire. If it is instead a mechanism to ensure better end joining via microhomology, this changes the paradigm.
* accurately measuring diversity in individual immune repertoires
* comparing repertoires between individuals 

## What are the risks and the payoffs?

There is a risk that we cannot obtain direct evidence of microhomology guiding V(D)J recombination events using repertoire sequencing data.
There is also a risk of a mixed signal, where microhomology partially guides V(D)J recombination events but is not the primary determinant.
This risk is mitigated by preliminary results we have obtained, and as described above, because this relationship has not been fully defined in a human system, any result (either positive or negative) will be informative.

The payoffs are in the applications above.

## Are there any other categorically different approaches that could be applied here?

I think our approach is the most basic one.
We could perhaps only have a data-driven approach, with no modeling, but it is not categorically different from ours.

## Is it possible that better data would make this project irrelevant?

We will never be able to manipulate the V(D)J recombination mechanism in vivo in humans.
Immune receptor repertoire data provide snapshots of the outcome of V(D)J recombination.
No other data type can currently show the underlying intermediary steps.

## Is there pre-existing work/code that could be leveraged to at least get a first-pass answer for the underlying scientific question? To validate the approach?

Previous work from Michael Lieber's group has suggested that microhomology can constrain junctional diversity during V(D)J recombination in-vitro (Gerstein and Lieber, Nature 1993), especially when TdT is not present at high levels (Gauss and Lieber, Mol. Cell. Biol. 1996).
While this work does not provide code or data that could be leveraged for our project, it does provide biological support for our hypothesis.

* Maggie's mechanistic trimming model work [here](https://github.com/magdalenarussell/mechanistic-trimming) 
* Maggie's work in the TCR-GWAS paper [here](https://github.com/phbradley/tcr-gwas)

Maggie's trimming model has identified some important sequence features that influence trimming probabilities; we can explore the relationship between these sequence features and microhomology.
Once we have a model, we can look for differences in microhomology between individuals and in the context of SNPs identified in the TCR-GWAS work. 

## What data will you use? Are there appropriate hold-out sets?

For model training and exploratory analysis:
* We will use simulated IGoR sequences
* We will use the Emerson data set as a very large, "typical blood repertoire" data set

For model validation/testing:
* We will use other "typical blood repertoire" TCR data sets, such as Paul Thomas and Aisha Souqette's data set
* We will use preterm and cord blood repertoires to limit N insertions
* We will use thymic TCR data from young children, for a bigger amount of non-productive sequences
* We will use some Omenn's Syndrome repertoires, where they have oligoclonal TCR repertoires (defective V(D)J recombination, will keep some of the model's "moving parts" stable)
* We will test on some BCR data to verify if the microhomology model applies across all immune receptors
* We will use the GWAS-Emerson set, to look for genomic variation associations

## Sketch the approach, broken down into steps, with expected amounts of time and intermediate steps for each.

### Preliminary signal check (1 months)--Assya:

Make sure there is signal in the data. 
* See if we see a signal for less trimming/N insertions in junctions with microhomology
* See if there is a bias in gene usage vs microhomology

TODO: add details of Assya's work thus far

See [here](https://github.com/magdalenarussell/mh-tex) for details.

#### Test

Do we see signal?

#### How do we decide to move onto the next stage?

If there is signal in the data, sounds reasonable to continue.

### Level 1 model (1 month):

Design a "mechanistic trimming" model to define the relationship between trimming and microhomology.
We will only model trimming at this stage, not N-insertions.
Ideally, this model will be structured in a way that will allow us to explore the relationship between sequence-level features, MH, and trimming probabilities.

We will restrict the training data set to include only non-productive sequences containing zero N-insertions.

#### Test

Does this model fit the no-N-insertion data better than Maggie's current "mechanistic trimming" model?

#### How do we decide to move onto the next stage?

Depending on the test result, we can move onto modeling sequences with N insertions.


### Level 2 model (2-3 months):

Design a "mechanistic" junctional processing model to define the relationship between trimming, N-insertion, and microhomology.
Here, we will build upon our model from Level 1 to add the N-insertion component.
Ideally, this model will be structured in a way that will allow us to:

- continue exploring the mechanistic relationship between trimming and MH
- explore the mechanistic relationship between N-insertion/other NHEJ processes and MH 
- explore whether the amount of MH influences the order in which nucleases/polymerases/ligases/etc. are recruited to the junction
- define the optimal amount of MH between two sequences during NHEJ 

We will continue restricting the training data set to include only non-productive sequences, but will now include all sequences, regardless of their extent of N-insertion.

#### Test

Does this model fit the no-N-insertion data better than the current existing model and the Level 1 model? What about the N-insertion data?

#### How do we decide to move onto the next stage?

Depending on the test result, we can move onto validating each model using various testing data sets.

### Level 3 (2 months):

Validate the models with testing data sets and write the manuscript.
