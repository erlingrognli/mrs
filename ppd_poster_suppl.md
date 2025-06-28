Supplementary material for the poster: “EOS and OCD - different brain
signatures are noted though overlapping features exist”
================
Nor C. Torp, Erling W. Rognli, Stener Nerland, Runar E. Smelror, Tord
Ivarsson & Ingrid Agartz
2025-06-28

### Posterior Predictive Distributions

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

![](C:/Users/erlrog/Documents/mrs/plots/area/ppd_area.png)

![](C:/Users/erlrog/Documents/mrs/plots/volume/ppd_volume.png)

[Return to main page](http://github.com/erlingrognli/mrs)
