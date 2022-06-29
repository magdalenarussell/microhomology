# Project overview

## What are you trying to do? Articulate your objectives using words that would be familiar to someone who had taken an undergraduate class in the topic.

Adaptive immune receptors are created via a stochastic process called V(D)J recombination. We are trying to verify whether microhomology guides V(D)J recombination events such as trimming and N nucleotide insertions.

## What are you curious about? What will keep pushing you forward when things get tough?

We're curious about how V(D)J recombination works!
Specifically, we're interested in whether microhomology is important for V(D)J recombination and if so, how many nucleotides of microhomology are optimal for the process.
We are also interested in the relationship between microhomology and trimming/N-insertions as well as the order in which they happen.
Since current V(D)J recombination models are not aware of microhomology, we want to incorporate microhomology information to improve model accuracy.

## How is it done today, and what are the limits of current practice?

1. The current V(D)J recombination models are not aware of microhomology. We expect that incorporating microhomology into model formulation will lead to a better fit to the data.
2. Currently, trimming and N insertion are assumed to be independent. If microhomology has an influence on these V(D)J recombination events, it will change the way these events are thought of.
3. V/D/J gene usage bias has been observed and reported (and learned by the current models described in point 1). We suspect that V/D/J gene usage is at least in part influenced by microhomology between the two joining segments. Therefore, a V/D/J bias might not need to be learned.

## What's new in your approach and why do you think it will be successful?

We will:
1. Formulate the microhomology-dependent V(D)J recombination model as a probabilistic model
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
However, securing collaborations for in vitro validation could be an option.
There is also a risk of a mixed signal, where microhomology partially guides V(D)J recombination events but is not the primary determinant.
This risk is mitigated by preliminary results we have obtained.

The payoffs are in the applications above.

## Are there any other categorically different approaches that could be applied here?

I think our approach is the most basic one.
We could perhaps only have a data-driven approach, with no modelling, but it is not categorically different from ours.

## Is it possible that better data would make this project irrelevant?

We will never be able to manipulate the V(D)J recombination mechanism in vivo in humans.
Immune receptor repertoire data provide snapshots of the outcome of V(D)J recombination.
No other data type can currently show the underlying intermediary steps.

## Is there pre-existing work/code that could be leveraged to at least get a first-pass answer for the underlying scientific question? To validate the approach?

Previous work from Michael Lieber's group has suggested that microhomology can constrain junctional diversity during V(D)J recombination (Gerstein and Lieber, Nature 1993), especially when TdT is not present at high levels (Gauss and Lieber, Mol. Cell. Biol. 1996).
While this work does not provide code or data that could be leveraged for our project, it does provide biological support for our hypothesis.

* Maggie's mechanistic trimming model work: https://github.com/magdalenarussell/mechanistic-trimming 
* Maggie's work in the TCR-GWAS paper: https://github.com/phbradley/tcr-gwas 

Maggie's trimming model has identified some important sequence features that influence trimming probabilities; we can explore the relationship between these sequence features and microhomology.
Once we have a model, we can look for differences in microhomology between individuals and in the context of SNPs identified in the TCR-GWAS work. 

## What data will you use? Are there appropriate hold-out sets?

%MR I think it would be good to clarify which data we want to train our model with (and use for exploratory analysis) and which we want to hold-out
%MR I have proposed a split here, but of course feel free to change whatever!
For model training and exploratory analysis:
* We will use simulated IGoR sequences.
* We will use preterm and cord blood repertoires to limit N insertions.
* We will use thymic TCR data from young children, for a bigger amount of non-productive sequences.

For model validation/testing:
* We will use some Omenn's Syndrome repertoires, where they have oligoclonal TCR repertoires (defective V(D)J recombination, will keep some of the model's "moving parts" stable)
* We will test on some BCR data to verify if the microhomology model applies across all immune receptors.
* We will use the GWAS-Emerson set, to validate in a large cohort and look for genomic variation associations

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
Identify relationship between trimming and microhomology.

#### Test

Does this model fit the no-N-insertion data better than the current existing model?

#### How do we decide to move onto the next stage?

Depending on the test result, we can move onto modeling sequences with N insertions.

### Level 2 model (2 months):

Add the N insertions component.
Identify relationship between trimming/N-insertion and microhomology.

#### Test

Does this model fit the no-N-insertion data better than the current existing model and the Level 1 model? What about the N-insertion data?

#### How do we decide to move onto the next stage?

Depending on the test result, we can move onto validating the model and improving speed.

### Poke at the model (2 months):

Look for someone who is willing to make this model fast.
Poke at the model.
Test on different held-out datasets.
Write manuscript draft.

#### Test
At this stage, if some things need to be verified in vitro, make a concrete plan of what needs to be verified and contact relevant teams for collaboration.
%MR Do we want to discuss how we will first explore the relationship between MH and trims/inserts and then explore how MH may influence gene usage?
%MR Or are you thinking it would be better to incorporate gene usage into the model from the start?
Train the model on equal amounts of each V/D/J pairing to make it not learn any V/D/J gene usage bias and how it performs
