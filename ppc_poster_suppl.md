Supplementary material for the poster: “EOS and OCD - different brain
signatures are noted though overlapping features exist”
================
Nor C. Torp, Erling W. Rognli, Stener Nerland, Runar E. Smelror, Tord
Ivarsson & Ingrid Agartz
2025-06-28

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
of the cortical area.
![](C:/Users/erlrog/Documents/mrs/plots/area/violin_ppc_area.png)

And similarly, we make the same type of plot for subcortical volumes:

![](C:/Users/erlrog/Documents/mrs/plots/volume/violin_ppc_volume.png)

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

![](C:/Users/erlrog/Documents/mrs/plots/area/loo_pit_area.png)

And for the model for subcortical volumes:

![](C:/Users/erlrog/Documents/mrs/plots/volume/loo_pit_volume.png)

Both also indicating adequate model fit.

[Return to main page](http://github.com/erlingrognli/mrs)
