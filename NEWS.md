
# spacemap 0.50.0

04/13/17

### Major changes

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
### Minor Changes

- Updated parameter names to have consistent camel case style with clearer names. 
