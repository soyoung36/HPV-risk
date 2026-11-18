# Third-party code

The LHS/PRCC routines were obtained from the University of Michigan Kirschner
laboratory: http://malthus.micro.med.umich.edu/lab/usadata/.

`LHS_Call.m` retains the original sampler and its source comments.
`PRCC_sen.m` was adapted for a single response and this model's parameter
selection, using MATLAB partial Spearman correlation. Original PRCC code
credits Simeone Marino (2007), with modifications by Marissa Renardy (2020).

Marino S, Hogue IB, Ray CJ, Kirschner DE (2008). A methodology for performing
global uncertainty and sensitivity analysis in systems biology.
*Journal of Theoretical Biology* 254, 178-196.
https://doi.org/10.1016/j.jtbi.2008.04.011

The sampling method is described in Stein M (1987), Large Sample Properties
of Simulations Using Latin Hypercube Sampling, *Technometrics* 29, 143-151.
Only the uniform sampling branch is used in this study.
