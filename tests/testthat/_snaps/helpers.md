# which_omics works

    Code
      which_omics()
    Message
      "rna-seq", also as "rnaseq", "rna", "gene counts", "transcriptomics"
          is modeled by default as count data with a negative-binomial marginal.
      
      "mirna-seq", also as "mirnaseq", "microrna-seq", "micrornaseq", "mirna",
          "microrna"
          is modeled by default as count data with a negative-binomial marginal.
      
      "proteomics", also as "protein fragments", "protein counts"
          is modeled by default as count data with a negative-binomial marginal.
      
      "metabolomics", also as "lc-ms", "gc-ms", "ms"
          is modeled by default as positive continuous data with a log-normal
          marginal.
      
      Anything else is modeled with the empirical marginal.

# which_marginals works

    Code
      which_marginals()
    Message
      "e" or "empirical" for using the empirical marginal distribution;
      "n" or "normal" for using the normal marginal distribution;
      "ln" or "lognormal" for using the log-normal marginal distribution;
      "nb" or "negative binomial" for using the log-normal marginal distribution.

