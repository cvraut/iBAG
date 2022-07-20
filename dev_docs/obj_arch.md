# Object Architectures
---
this is a help doc created to organize key architecture designs

*NOTE*: using an object oriented approach so that we can construct new versions
of iBAG like surviBAG or piBAG using inheritance & we can establish a consistent
 framework.

----------
## objects.R
----------
iBAG_data
  - attributes:
    - supplied (default):
      - mrna (demo_mrna)
      - outcome (demo_outcome)
      - data_1+ (data_1=demo_cnv)
      - DEBUG (FALSE)
      - one_val_per_gene (TRUE)
    - inferred:
      - n_genes
      - n_patients
      - n_data
      - genes
      - patients
      - data_names
  - methods:
    - update_data
    - validate
    - getters
      - patients
      - genes
      - n_patients
      - n_genes
      - mrna
      - outcome
      - one_val_per_gene
      - data
      - data_names

iBAG_results
  - attributes:
    - supplied (default):
      - X (NULL)
      - Y (NULL)
      - SS (list(SST = NULL, SSE = NULL))
      - beta_mean (NULL)
      - beta_incl_prob (NULL)
  - methods:
    - getters
      - X
      - Y
      - SS
      - beta_mean
      - beta_incl_prob
    - setters
      - X
      - Y
      - SS
      - beta_mean
      - beta_incl_prob
