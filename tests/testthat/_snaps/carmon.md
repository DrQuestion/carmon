# verbosity 2 and 0 work

    Code
      c_obj <- carmon(multi_omics_micro, net_method = "correlation", cor_cutoff = 0.7,
        verbose = 2, plot = FALSE)
    Message
      Checking sample-matching, formatting the data,
          and other formalities....
      Done!
      
      ****************Beginning copulization****************
      Copulizing layer 1 of 2 (rnaseq)
      Using negative binomial marginals....
      Rnaseq layer copulized!
      Copulizing layer 2 of 2 (metabolomics)
      Using lognormal marginals....
      Metabolomics layer copulized!
      ****************Copulization complete*****************
      
      ***********Beginning network reconstruction***********
      Reconstructing network with Pearson's correlation....
      Correlation cutoff is 0.7
      Network reconstructed!
      ***********Network reconstruction complete************
      
      **************Beginning network analysis**************
      Computing centrality measure 1 of 4 (Degree)
      No node was found to be central for Degree.
      Computing centrality measure 2 of 4 (Betweenness) 
      No node was found to be central for Betweenness.
      Computing centrality measure 3 of 4 (Closeness) 
      No node was found to be central for Closeness.
      Computing centrality measure 4 of 4 (Eigenvector) 
      No node was found to be central for Eigenvector Centrality.
    Condition
      Warning in `compute_centrality()`:
      No node was found to be central. No centrality report
          will be generated.
    Message
      **************Network analysis complete***************

# verbosity 1 works

    Code
      c_obj <- carmon(multi_omics_micro, net_method = "correlation", cor_cutoff = 0.7,
        verbose = 1, plot = FALSE)
    Message
      Checking sample-matching, formatting the data,
          and other formalities....
      Done!
      
      ****************Beginning copulization****************
      Copulizing layer 1 of 2 (rnaseq)
      Copulizing layer 2 of 2 (metabolomics)
      ****************Copulization complete*****************
      
      ***********Beginning network reconstruction***********
      Reconstructing network with Pearson's correlation....
      Correlation cutoff is 0.7
      Network reconstructed!
      ***********Network reconstruction complete************
      
      **************Beginning network analysis**************
      Computing centrality measures....
      No node was found to be central for Degree.
      No node was found to be central for Betweenness.
      No node was found to be central for Closeness.
      No node was found to be central for Eigenvector Centrality.
      Centrality measures computed!
    Condition
      Warning in `compute_centrality()`:
      No node was found to be central. No centrality report
          will be generated.
    Message
      **************Network analysis complete***************

# print.carmon without centrality works

    Code
      print(c_obj)
    Output
      Carmon network estimated with correlation.
      
      The call was:
      carmon(layers = multi_omics, net_method = "correlation", cor_quant = 0.05, 
          analyse = FALSE, plot = FALSE, verbose = 0)
      
      ******************************************************
      
      The network is made of 2 omics layers: rnaseq and metabolomics.
      It has a total of 238 nodes.
      The rnaseq layer has 162 nodes; 
          its chosen marginal distribution is the negative binomial.
      The metabolomics layer has 76 nodes; 
          its chosen marginal distribution is the lognormal.
      
      ******************************************************
      
      Find central nodes with:
      c_obj <- compute_centrality(c_obj)
      Then print a report of the central nodes with:
      centrality_report(c_obj)
      
      Plot a comparison of the central nodes with:
      plot_report(c_obj)
      
      Plot the carmon network with:
      plot(c_obj)

# print.carmon with centrality works

    Code
      print(c_obj)
    Output
      Carmon network estimated with correlation.
      
      The call was:
      carmon(layers = multi_omics, net_method = "correlation", cor_quant = 0.05, 
          plot = FALSE, verbose = 0)
      
      ******************************************************
      
      The network is made of 2 omics layers: rnaseq and metabolomics.
      It has a total of 238 nodes.
      The rnaseq layer has 162 nodes; 
          its chosen marginal distribution is the negative binomial.
      The metabolomics layer has 76 nodes; 
          its chosen marginal distribution is the lognormal.
      
      ******************************************************
      
      Print a report of the central nodes with:
      centrality_report(c_obj)
      
      Plot a comparison of the central nodes with:
      plot_report(c_obj)
      
      Plot the carmon network with:
      plot(c_obj)

