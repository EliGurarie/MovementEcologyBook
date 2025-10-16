#' ---
#' title: "Segmentation Part III: Hidden Markov Models"
#' author: "Elie Gurarie"
#' subtitle: "EFB 796: Techniques in Movement Ecology"
#' date: "October 2, 2025"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#' ---
#' 
## ----setup, include=FALSE----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, 
                      warning = FALSE, cache = TRUE)

#' 
#' # Goals
#' 
#' > To estimate a Multi-State Random Walk for some elk data, using the `momentuhmm` package. 
#' 
#' Useful packages:
#' 
## ----------------------------------------------------------------------------------------------------
require(plyr)
require(mapview)
require(sf)
require(ggplot2)

#' 
#' 
#' And load the *regularized* elk data here: [elk_regdata.rda](elk_regdata.rda). 
#' 
## ----echo = -1---------------------------------------------------------------------------------------
load("elk_regdata.rda")
head(elk_regdata)

#' 
#' 
#' 
#' # Select Individuals (Non-Residential)
#' 
#' We can one of the visualization methods we learned to examine the tracks for non-residential (migratory) behavior.
#' 
## ----fig.height=4, fig.width=6-----------------------------------------------------------------------
ggplot(data=elk_regdata, aes(x=lon, y=lat)) +
  geom_path(size=0.5) +
  theme_classic() +
  facet_wrap(~id, scale="free", ncol=4)

#' 
#' This subset of animals has fairly clear migratory movements.  We'll see if we can parse those movements even more finely. 
#' 
#' 
#' NOTE: this version of the elk data has been "regularized" at 2 hour intervals.  Here's one way to see that:
#' 
## ----------------------------------------------------------------------------------------------------
dlply(elk_regdata, "id", 
      function(df) summary(diff(df$date)))

#' 
#' 
#' # Hidden Markov Model
#' 
#' Hidden Markov Models (HMMs) can be used for time series data to describe 2 processes, one observed (e.g., the data produced by telemetry, such as locations and movement metrics) and one "hidden", representing the "hidden" ecological behaviors driving the observed data.
#' 
#' The Markov chain portion of HMM's state that the probability of being in a particular behavior state at time $t+1$ is **only** dependent on the current state at time $t$. Once these probabilities are described, a transition probability matrix, or the probability of moving between behavior states, can be described. HMM's can have any number of hidden behavioral states, granted that they are biologically meaningful (model selection can also assist with determining the "optimal" number of states).
#' 
#' Importantly, HMMs are performed in discrete time and assume that observations are "conditionally independent". These can be difficult for movement data derived from tags, which can be prone to irregularly spaced observations and errors.
#' 
#' ## Load the data
#' 
#' Now that our data is regular space/time, we can use the methods we learned in lab 2 to make our data spatial, convert the coordinates to a projected coordinate system, and save the X/Y coordinates as columns in our data.
#' 
## ----------------------------------------------------------------------------------------------------
library(sf)

#' 
#' Prepare 
#' 
## ----------------------------------------------------------------------------------------------------
myelk_utm <- elk_regdata |> subset(id == "YL96") |> 
  st_as_sf(coords = c("lon","lat"), crs = 4326) |>
  st_transform(32611)

#' 
#' We will need the coordinates!  But the data should be a data frame
#' 
## ----------------------------------------------------------------------------------------------------
myelk_utm <- cbind(myelk_utm, st_coordinates(myelk_utm)) |> data.frame()
head(myelk_utm)

#' 
#' ## Fit HMM: Two States
#' 
#' Now that our data is regularized and we have our movement metrics annotated to our data, we can fit our HMM with the `momentuHMM` package, which is very flexible.  You can use  user-specified parameter distributions, options for hierarchical and multivariate models, biased correlated random walks and more. As usual, read the [vignette](https://cran.r-project.org/web/packages/momentuHMM/vignettes/momentuHMM.pdf).
#' 
## ----------------------------------------------------------------------------------------------------
library(momentuHMM)

#' 
#' We first prepare our data using the `prepData` function to grab our X/Y coordinates.
#' 
## ----------------------------------------------------------------------------------------------------
myelk_hmm_prep <- prepData(myelk_utm |> data.frame(), 
                           type = "UTM", coordNames = c("X","Y")) 

head(myelk_hmm_prep)

#' 
#' Note  that once again - this package made a data frame that computes all the things we computed, and `ltraj` computes. 
#' 
#' It is important that step lengths be greater than 0:
#' 
## ----------------------------------------------------------------------------------------------------
table(myelk_hmm_prep$step == 0)

#' 
#' MomentuHMM fits HMM's within a Bayesian framework, where prior distributions can be specified for each state-specific parameter to help distinguish the "true" behavioral state when combined with observed distributions. 
#' 
#' We will start with a two-state model with state-specific and distribution specific parameters for each variable that will be used to distinguish our behavioral states, namely the step length and relative turning angles.
#' 
#' We will make our priors somewhat informative, assuming that one state will be defined by smaller step lengths (smaller mean and sd) and bigger turning angles (smaller concentration parameter), with the opposite for the second state. 
#' 
## ----------------------------------------------------------------------------------------------------
stepMean0 <- c(m1 = 100, m2 = 4000)
stepSD0 <- c(sd1 = 50, sd2 = 1000)
angleCon0 <- c(rho1  = 0.1, rho2 = 0.8)

#' 
#' We will now assign names to our states, with the slower, more tortuous state being defined as a "resident" state and the faster, straighter state as a "transit" state.
#' 
## ----------------------------------------------------------------------------------------------------
stateNames <- c("resident","transit")

#' 
#' We now pick distributions for our variables of interest (step length and turning angle), using a Gamma distribution for the step lengths (`?dgamma`) and a Wrapped Cauchy distribution for the turning angles (`dwrpcauchy`).
#' 
#' We store our priors defined in the code chunk above into a list object, "Par0".
#' 
## ----------------------------------------------------------------------------------------------------
dist <- list(step = "gamma", angle = "wrpcauchy")
Par0 <- list(step=c(stepMean0, stepSD0), angle = c(angleCon0))

#' 
#' We can now fit our HMM using the `fitHMM` function with our prepped data and objects from above. We also specify the number of states to identify (`nbStates`).
#' 
## ----------------------------------------------------------------------------------------------------
myelk_hmm_fit <- fitHMM(data = myelk_hmm_prep, nbStates = 2, dist = dist, 
                      Par0 = Par0, stateNames = stateNames)

print(myelk_hmm_fit)

#' 
#' After fitting our model, we can use the Viterbi formula with the `viterbi` function to decode the transition probabilities for the most likely states at each time point along our movement trajectory.
#' 
## ----------------------------------------------------------------------------------------------------
hmm_states <- viterbi(myelk_hmm_fit)
str(hmm_states)

#' 
#' We can add our predicted states as a column for each observation to our data:
#' 
## ----------------------------------------------------------------------------------------------------
myelk_utm$state <- hmm_states

#' 
#' ### Plot 2-state predictions
#' 
#' The default plotting of `momentuhmm` is quite useful: 
#' 
## ----DefaultTwoStatePlot, echo  = 2, eval = 1--------------------------------------------------------
plot(myelk_hmm_fit, ask = FALSE)
plot(myelk_hmm_fit)

#' 
#' 
#' We can then plot the results, using Base R to make multidimensional plots of the trajectory with the annotated states:
#' 
## ----TwoState_ScanTtrack-----------------------------------------------------------------------------
layout(cbind(c(1,1),2:3))
par(bty = "l", mar = c(2,2,2,2))

with(myelk_utm, {
  plot(X, Y, asp =1, col = c("orange","blue")[state], pch = 19, cex = 0.7)
  segments(X[-length(X)], Y[-length(Y)], 
           X[-1], Y[-1], col = c("orange","blue")[state[-length(state)]])
  plot(date, X, col = c("orange","blue")[state], pch = 19, cex = 0.7)
  segments(date[-length(X)], Y[-length(Y)], 
           date[-1], Y[-1], col = c("orange","blue")[state[-length(state)]])
  plot(date, Y, col = c("orange","blue")[state], pch = 19, cex = 0.7)
  segments(date[-length(X)], Y[-length(Y)], 
           date[-1], Y[-1], col = c("orange","blue")[state[-length(state)]])
})

#' 
#' We can also convert our data to an sf structure and use the methods we learned in lab 1 to plot the trajectory with each location colored by the predicted state with the mapview package.
#' 
## ----------------------------------------------------------------------------------------------------
library(mapview)

#' 
#' 
## ----mapviewTwoStateElk, cache = FALSE---------------------------------------------------------------
myelk_sf <- myelk_utm |>
  st_as_sf(coords=c("X","Y"), crs= 32611) |>
  st_transform(4326) |>
  mutate(state = as.character(state))

myelk_track <- myelk_sf |>
  dplyr::summarize(do_union=FALSE) |> 
  st_cast("LINESTRING")

require(mapview)
mapview(myelk_track, color="darkgrey") +
  mapview(myelk_sf, zcol="state", 
          col.regions=c("orange","blue"))

#' 
#' We can see from the plot that our movement "blobs" still have a mixture of both colors. We can try fitting a multivariate HMM, with latitude as a covariate, to assist with further separating our two states along our movement track.
#' 
#' This may suggest that there are actually 3 behavioral states here: one that is fast and straight (blue state), one that is slower and tortuous (orange state), and a third that is intermediate, where the animal is having directed, faster movements within its residential patch.
#' 
#' 
#' ## Fit HMM: 3-State Model
#' 
#' We can fit a 3 state model by simply adding an additional state-specific prior for each movement variable.
#' 
## ----Priors_ThreeState-------------------------------------------------------------------------------
stepMean0 <- c(m1 = 50, m2 = 2000, m3 = 200)
stepSD0 <- c(sd1 = 50, sd2 = 1000, sd3 = 100)
angleCon0 <- c(rho1  = 0.1, rho2 = 0.8, rho3 = 0.2)

#' 
#' We will label this third state as an additional "faster" residential state.
#' 
## ----------------------------------------------------------------------------------------------------
stateNames <- c("resident-slow","transit", "resident-faster")

#' 
## ----------------------------------------------------------------------------------------------------
dist <- list(step = "gamma", angle = "wrpcauchy")
Par0 <- list(step=c(stepMean0, stepSD0), angle = c(angleCon0))

#' 
#' We re-fit the model as before, but this time specifying 3 states.
#' 
## ----Fit_ThreeState----------------------------------------------------------------------------------
myelk_hmm_threestate <- fitHMM(data = myelk_hmm_prep, nbStates = 3, 
                         dist = dist, Par0 = Par0, 
                         stateNames = stateNames, 
                         formula = ~1)

#' 
## ----------------------------------------------------------------------------------------------------
myelk_hmm_threestate

#' 
#' 
#' We can decode our model results and bind the predicted states for each observation back to our data:
#' 
## ----------------------------------------------------------------------------------------------------
hmm_3states <- viterbi(myelk_hmm_threestate)
myelk_utm$state <- hmm_3states

#' 
#' ### Plot 3-sate predictions
#' 
#' Now if we plot the results with our new state in green, we can see that our model did a much better job at finding our transiting state and two residential states within our movement "blobs"!
#' 
## ----ThreeStateDefaultPlot, echo = 2, eval = 1-------------------------------------------------------
plot(myelk_hmm_fit, ask = FALSE)
plot(myelk_hmm_fit)

#' 
#' Here is a mapview: 
#' 
## ----ThreeStateMapview, cache = FALSE----------------------------------------------------------------
myelk_track <- myelk_sf |>
  dplyr::summarize(do_union=FALSE) |> 
  st_cast("LINESTRING")
myelk_sf <- myelk_sf |> mutate(state = hmm_3states)
mapview(myelk_track, color="darkgrey") +
  mapview(myelk_sf, zcol="state", col.regions=c("orange","blue", "green"))

#' 
#' 
