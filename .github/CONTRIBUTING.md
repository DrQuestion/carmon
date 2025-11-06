# Contributing to carmon

This outlines how to propose a change to carmon.
For a detailed discussion on contributing to this and other tidyverse packages, please see the [development contributing guide](https://rstd.io/tidy-contrib) and our [code review principles](https://code-review.tidyverse.org/).

## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file. 
This generally means you'll need to edit [roxygen2 comments](https://roxygen2.r-lib.org/articles/roxygen2.html) in an `.R`, not a `.Rd` file. 
You can find the `.R` file that generates the `.Rd` by reading the comment in the first line.

## Bigger changes

If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it’s needed. 
If you’ve found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).
See our guide on [how to create a great issue](https://code-review.tidyverse.org/issues/) for more advice.

### Adding new omics types and marginal distributions

The carmon package is built to grow. We primarily plan in the future to expand 
the range of omics types carmon is tailored for. We also plan to expand the 
range of marginal probability distributions available. If you want to contribute
to this growth, please respect the inner function structure and environment of 
the package. We describe main guidelines below for a few different scenarios. 
Please, respect what is explained below, at the same time complying with the 
remaining rules for contribution. In any case, feel free to issue a feature 
request on the GitHub page of the package if you believe we should add support 
for new omics types or marginals.

#### New omics type, marginal distribution already implemented

If the marginal distribution you would want to assign as a default to your new
omics type is already available (check which ones are with `which_marginals()`),
these are the passages to follow to include a new omics type:

*   Add the omics type and a series of synonyms you believe appropriate to the 
    function `omics2marginal()` in the `R/copulize.R` file. For example, you 
    can see how "rna-seq" and "gene count" both link to the Negative Binomial 
    marginal distribution.

*   Modify the user helping function `which_omics()` and the input check 
    function `check_omics()` in the `R/helpers.R` file accordingly, adding your 
    new omics type, the accepted synonyms, and, in the case of `which_omics()`,
    the default marginal distribution you decided to assign to it. Please, 
    respect the syntax and structure (indentation, order of words, etc.) of both
    functions **and** of the message returned by `which_omics()` to mirror the 
    behaviour of the omics types that are already present in those.
    
#### New marginal distribution, no new omics type

If you want to expand the library of marginal distributions of the package, say 
you want to add the **F**oo **B**ar marginal distribution:

*   Create the function **fb**`_m_copulizer()` (notice the initials of the 
    distribution) in the `R/copulize.R` file. The function should estimate the 
    parameters of the marginal distribution from the layer it is given as an 
    input, using the best estimation strategy to your knowledge. 
     * If possible, it should include support for a model design matrix and 
       design formula to include covariates in the estimation of such 
       parameters. (*Not yet implemented as of version 0.99.0*)
     * It should also include strategies to perform inversion of non invertible 
       distributions. Implementing at least the median method is mandatory. 
     * Make sure to respect all the dimensionalities expected by the functions 
       you will use internally to the new function you are implementing. A good 
       way to check for it is to make sure that every variable in the layer you 
       obtain as an output is approximately distributed according to a standard
       normal distribution.
     * Make sure to appropriately deal with occasional infinite values with the 
       function `hackInf()`.
     * Follow the already existing functions as examples. The `nb_m_copulizer()`
       function as an example of a non-invertible marginal, and the 
       `ln_m_copulizer()` as an example of an invertible one.

*   Add an `else if` statement to the function `copulize()` in the file 
    `R/copulize.R`, where indicated, inside the for loop where each layer 
    is copulized, where your `fb_m_copulizer()` function will be used.
    *Follow the same syntax already in use for the other marginals*,
    and everything should work fine. In particular, make sure to include both 
    the initials and the full marginal distribution name as in the provided 
    examples. In this case, it would be with `"fb`" and `"foo bar"`, as follows:
     * `else if (tolower(marginals[l] %in% c("fb", "foo bar")) { ... }`
     * Make sure to overwrite the `marginals[l]` element with the `"foo bar"`
       full name.
       
*   Modify the `which_marginals()` function in the `R/helpers.R` file
    to include in the message information about your now marginal, following the
    same syntax and indentation of the ones that are present already.
    
#### New omics type, new marginal distribution

If you are adding a new omics type which follows a marginal distribution that is
not implemented yet in the `carmon` package, please follow all the steps above 
apart from the very first step described, the one modifying the 
`omics2marginals()` in the `R/copulize.R` file. Instead, modify it as follows:

*   To allow the new marginal distribution to be used as a default for the new 
    omics type, modify the `omics2marginals()` function by adding an `else if` 
    statement **above** the final else statement that assigns the empirical 
    marginal distribution. Also in this case, include the omics type and a 
    series of synonyms you believe appropriate, and assign the new marginal 
    distribution to `marginals[l]` object with `marginals[l] <- "fb"`. Again, 
    *following the same syntax already available should allow this change to be*  
    *accomplished seamlessly*.

### Pull request process

*   Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis::create_from_github("DrQuestion/carmon", fork = TRUE)`.

*   Install all development dependencies with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
    If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a Git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

*  For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just below the first header). Follow the style described in <https://style.tidyverse.org/news.html>.

### Code style

*   New code should follow the tidyverse [style guide](https://style.tidyverse.org). 
    You can use [Air](https://posit-dev.github.io/air/) to apply this style, but please don't restyle code that has nothing to do with your PR.  

*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation.  

*  We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. 
   Contributions with test cases included are easier to accept.  

## Code of Conduct

Please note that the carmon project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
