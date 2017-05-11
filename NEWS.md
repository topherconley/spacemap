# spacemap 0.53.0

04/13/17

Passing R CMD check with no errors. Small warning on LICENSE in DESCRIPTION.

### Minor changes

- Updated .Rd files to pass R CMD check.
- Pulled in Jie's edits of vignettes.
- Resolved NAMESPACE issues by using Roxygen2.
- Adjusted Depends: and Imports:  fields in DESCRIPTION 
according to best practices.

# spacemap 0.50.0

04/13/17

### Minor changes

- Input parameters and output for functions `cvVote`, `bootEnsemble`, 
`bootVote` has been restructured to be more standardized with previous
parameters given defaults and only able to change through ... construct. 
- Edited vignette basics.Rmd to be more clear. 
- Created vignettes tuning.Rmd and ensemble.Rmd. 
Removed old tuning_sim1.html from docs/reference. 
- More accurate documentation for model fitting, model selection functions. 
Deleted documentation that should not be exposed to user. 
- For function `space.joint` and `spacemap`, if sig and rho are non-null 
parameters, they are still updated. 
- Updated parameter names to have consistent camel case style with clearer names. 
