# Bachelor thesis
This Github is part of the bachelor thesis by Yannick van Etten. 
Several methods are evaluated to reduce the estimation errors of the weights in forecast combinations. Radchenko et al. (2023) evaluate a trimming method in detail. 
In the thesis, this will be replicated. This is a way to diminish the impact of estimation errors. Furthermore, the thesis tries to find other ways to reduce the presence of 
estimation errors. The impact of using different covariance estimation methods is explored. The thesis can be read [here](https://github.com/YannickvanEtten/bachelor_thesis/blob/7fb650b04c39e92f6c4872ae7b096a1e0a80203f/Bachelor_thesis_Yannick_van_Etten_final.pdf).

This first Notebook, Simulation_study_trimming_methods, focuses on Figure 6 from Radchenko et al. (2023), where the mean squared forecasting error (MSFE) 
of a forecast combination is plotted conditional to a trimming threshold. In addition, their paper looks at the impact of fixed weights on the MSFE. Lastly, 
this Notebook recreates Figure 7 from Radchenko et al. (2023). In this plot, the MSFE of a forecast combination is plotted conditional to a trimming threshold, 
similar to Figure 6. However, this time this is done for different AR(1) time series, changing the autoregressive parameter $\phi_1$. The Notebook can be viewed [here](https://github.com/YannickvanEtten/bachelor_thesis/blob/1c83e666738d030138d2244ad669f6fca952c85e/Simulation_study_covariance_estimation.ipynb).

The second Notebook, Simulation_study_covariance_estimation, focuses on the performance of different covariance estimation methods. Especially the linear shrinkage method and the factor model method. These methods are evaluated on the hand of the MSFE of a forecast combination. The trimming method is applied to the weights obtained using 
these alternative covariance estimation methods. The Notebook can be viewed [here](https://github.com/YannickvanEtten/bachelor_thesis/blob/1c83e666738d030138d2244ad669f6fca952c85e/Simulation_study_trimming_method.ipynb).

Lastly, a file with R code is added to this page. This contains the data-driven threshold method code from Radchenko et al. (2023). In their paper, they make a mistake in calculating the p-values. This mistake is solved in this file.
