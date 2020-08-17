# Power/Sample Size App for AD Clinical Trials

# Study designs/app sections:

- MMRM (basic)
- MMRM with ADNI-based pilot estimate generator 

# Inputs (basic): (allow range slider for one?)

- n number of observations (per group)
- delta true difference in means
- ...

# Inputs (ADNI pilot generator): 

- inc/exc criteria
- outcome
- duration
- ...

# Outputs:

- plot with omitted input on y-axis and range input on the x-axis
- some text summary

## Two tabs or two separate apps
- Basic interface to longpower::mmpower.mmrm (assuming user already know all pilot parameters.
- Use ADNIMERGE::adnimerge to generate pilot estimates
  - Challenge: mismatch ADNI followup and new trial design followup schedule
    - Pilot models should use continuous-time where possible. Most natural way to extrapolate (with limitations, e.g. learning effects [clear caveat])
    - Random slope models?
    - gls correlation structure CAR; heterogeneous variance?
- Duration of followup?
  - The new trial design has to limited to available followup in ADNI?
  - Maybe warning where extrapolating off of ADNI support?
- At any rate, we should read out some summary of pilot model fit
- Inc/Exc criteria:
  - Age ranges, APOE4, AV45 range, ABETA range, Tau,  Tau/Abeta (Roche Elecsys platform)
  - Baseline Dx, MMSE range, CDR-G, CDR-SB, Logical Memory (LDELTOTAL)
  - Shiny double slider
- Covariate options:
  - Age, gender, education, eth, race, apoe4, AV45
- Outcome options:
  - FDG PET, AV45 PET, ADAS11, ADAS13, MMSE, CDRSB, mPACCtrailsB
- Automated model checking (enough observations, residuals)
- Imputation?
- Describe/present/plot ADNI sample used; and model fit
- ADNI-based estimates of attrition 
