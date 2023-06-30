# bachelor_thesis
This Github is part of the bachelor thesis by Yannick van Etten. 
Several methods are evaluated to reduce the estimation errors. In Radchenko et al. (2023) a trimming method is evaluated in detail. In the thesis this will be replicated. 

This first Notebook, Simulation_study_radchenko, focusses on the Figure 6 from Radchenko et al. (2023), where the mean squared forecasting error (MSFE) 
of a forecast combination is plotted conditional of a trimming threshold. In addition, their paper looks to the impact of fixed weights on the MSFE. Lastly, 
this Notebook recreates Figure 7 from Radchenko et al. (2023). In this plot the MSFE of a forecast combination is plotted conditional of a trimming threshold, 
similar to Figure 6. However, this time this is done for different AR(1) timeseries, changing the autoregresive paramater $\phi_1$. 

The second Notebook, Simulation_study_extension, focusses on the performance of different covariance estimation methdods. Specially, the linear shrinkage method and the 
factor model method. These methods are evaluated on the hand of the MSFE of a forecast combination. The trimming method is applied on the weights obtain by the usage of 
these alternative covariance estimation methods.
