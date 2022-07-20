# iBAG: integrative Bayesian Analysis of the Genome
An R package for the iBAG library

## About iBAG

The iBAG framework allows us to study the genomic effects from multiple types of genomic data in a hierarchical fashion.

![iBAG diagram](https://cvraut.github.io/iBAG_supplementary/images/iBAG_intro_pic.PNG "iBAG figure 1")

## Get started with iBAG in 2 lines:
```R
demo_data <- iBAG::demo_data()
iBAG_result <- iBAG::fit(demo_data)
```
```R
summary(iBAG_result)
>
summary(iBAG_result$mechmodel)
>
```