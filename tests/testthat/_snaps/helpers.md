# which_omics work

    Code
      which_omics()
    Output
      "rna-seq", also as "rnaseq", "gene counts", "transcriptomics"
         is modeled by default as count data with a negative-binomial marginal.
      
      "proteomics", also as "protein fragments", "protein counts"
         is modeled by default as count data with a negative-binomial marginal.
      
      "metabolomics", also as "lc-ms", "gc-ms", "ms"
         is modeled by default as positive continuous data with a log-normal marginal.
      
      Anything else is modeled with the empirical marginal.

# which_marginals work

    Code
      which_marginals()
    Output
      "e" or "empirical" for using the empirical marginal distribution;
      "n" or "normal" for using the normal marginal distribution;
      "ln" or "lognormal" for using the log-normal marginal distribution;
      "nb" or "negative binomial" for using the log-normal marginal distribution.

