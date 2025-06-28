Supplementary material for the poster: “EOS and OCD - different brain
signatures are noted though overlapping features exist”
================
Nor C. Torp, Erling W. Rognli, Stener Nerland, Runar E. Smelror, Tord
Ivarsson & Ingrid Agartz
2025-06-28

## Prior distributions in the model

## Posterior predictive checking

To evaluate model fit, we can use posterior predictive checking. If the
model fits the data well, it should be able to predict the outcome
variable from the predictor variables and the fitted model parameters
with reasonable accuracy. We can evaluate this by making repeated draws
of the predicted distribution of the outcome variable given the
predictor variables and the fitted model, and then plotting these
against the observed outcome variable that was used to fit the model in
the first place. Model fit can then be evaluated directly from the
similarity of the distributions.

Here we have made violin plots for this purpose, using the excellent
*bayesplot* package. The blue, shaded distributions are plotted from
4000 draws of the outcome variable predicted by the model, while the
dots and the violin is the observed distribution. We plot the various
brain regions separately. Note that the scale on the y-axis is the log
of the cortical area. ![](plots/area/violin_ppc_area.png)

And similarly, we make the same type of plot for subcortical volumes:

![](plots/volume/violin_ppc_volume.png)

As is evident from the plots, the model fits reasonably well.

Another form of posterior predictive checking available in *bayesplot*
is the loo-pit plot. This uses approximate leave-one-out
cross-validation (through the excellent *loo* package) and the
probability integral transform, calculating the cumulative probability
density of each observation when left out of fitting the model. With a
calibrated model, these follow a uniform distribution between 0 and 1
(see bayesplot documentation for references). The following plot
compares the distribution of these loo-pit values to repeated draws from
a uniform 0-1 distribution. Deviations from the uniform distribution can
be used to diagnose unmodelled under- or overdispersion, or bias.

The loo-pit plot for the model for cortical areas:

![](plots/area/loo_pit_area.png)

And for the model for subcortical volumes:

![](plots/volume/loo_pit_volume.png)

Both also indicating adequate model fit.

## Posterior Predictive Distributions

Looking at coefficient estimates in isolation is not always a good way
of understanding the implications of a fitted model. However, given any
set of predictor variables we choose, we can use the joint posterior
distribution to draw the implied distribution of the modelled variable.
This is called the posterior predictive distribution - the predicted
distribution of the outcome variable conditional on the fitted model and
the chosen set of predictor variable values. Plotting the posterior
predictive distributions is particularly illustrative.

Here we plot the posterior predictive distributions for EOS, OCD and
Controls separately. When simulating these distributions we have used
the proportion of females in the sample as the predictor value for
female gender, and the sample mean of age as the predictor value for
age. We have omitted the random intercepts per participants, as these
were based on the mean-centered values of intracranial volume. These
plots can hence be interpreted as the distribution of areas or volumes
for new controls, EOS-patients and OCD-patients, of indeterminate age
and gender.

As they are based on the posterior distribution, and not on point
estimates, these plots contain both the modelled variability (the
residual variance of the model) and the statistical uncertainty of all
the estimates.

![](plots/area/ppd_area.png)

![](plots/volume/ppd_volume.png)

[Return to main page](http://github.com/erlingrognli/mrs)
