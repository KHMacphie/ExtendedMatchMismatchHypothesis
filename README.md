# ExtendedMatchMismatchHypothesis
Data and code for "Timing isn’t everything: impacts of maximum abundance and duration of a seasonal resource on consumer fitness"

ACCESS INFORMATION
1. Licenses/restrictions placed on the data or code
NA
2. Data derived from other sources
NA
3. Recommended citation for this data/code archive
Will be updated after publication

DATA & CODE FILE OVERVIEW

This data repository consist of 5 data files, 12 code scripts, and this README document, with the following data and code filenames and variables

Data files and variables
1. Adults.csv: Adults breeding at each nestbox,
          site = site code,
          year = year,
          season = spring or winter catching (winter removed),
          box = nestbox number,
          ring = unique ID ring on bird,
          sex = sex of bird (M/F),
2. Bird_Phenology.csv: Each nest box in each year - breeding phenology and success,
          year = year,
          site = site code,
          box	= nextbox number, 
          species	= species breeding in box (coal tits removed),
          fki	= first known incubation date (ordinal 1=1st Jan),
          hd_1.45 = estimated hatch date,
          hatching_first_recorded	= date hatch date was first recorded,
          suc = how many successfully fledged,
3. Branch_Beating.csv: Caterpillar sampling data,
          site = site code,	
          year = year,
          date = ordinal date of each sample (1=1st Jan),
          tree = tree ID number,
          caterpillars = number of caterpillars recorded,	
          recorder = initials of sample recorder,
4. clutchswaps.csv: Nests included in clutch swaps (removed from analyses),
          site = site code,
          year = year,
          origin.nest	= nest where clutch originated,
          destination.nest = nest where clutch was reared,
5. Nestlings.csv: Individual nestling survival,
          site = site code,
          box	= nestbox number,
          year = year,	
          fledged = whether the nestling successfully fledged (0/1)
    
Code scripts and workflow

Run in order (OrganiseDataframes.R and ResultsFunctions.R are sourced in scripts)
1. OrganisingDataframes.R = Combining relevant data files and extracting required fields
2. EMHH_Model.R = RStan code for main Fitness Model of the EMMH
3. ResultsFunctions.R = Functions required for results scripts
4. FitnessModel_SummaryTable.R = Summary table of Fitness Model
5. FitnessFunction_Predictions.R = Predictions of fitness functions by caterpillar distribution
6. Slopes_FitnessFunctions_Plots.R = Plotting the effects of the caterpillar distribuiton on the fitness function
7. Lag_Predictions.R = Calculating and plotting predictions of lag, inc. Hatch Date Model
8. Selection_MeanFitness_Average.R = Estimates of average mean fitness components and selection on hatch date
9. Selection_MeanFitness_byLag.R = Estimates of varying mean fitness componenets and selection with lag
10. Selection_MeanFitness_Spatiotemporal.R = Site-year estimates of mean fitness and selection and their spatiotemporal decomposition
11. EMMH_ModelToCheckFit.R = RStan code for model to check data fit, predicting missing site-years
12. CheckFit_Simulations.R = Simulating under model to check fit to data

SOFTWARE VERSIONS

R Core Team. (2023). R v4.3.1: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.

Stan Development Team. (2022). RStan: the R interface to Stan. R package version 2.26.13. http://mc-stan.org/.

REFERENCES

NA
