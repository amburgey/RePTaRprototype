# RePTaR Prototype

We conceptualized, created, and trialed a prototype of a remote passive integrated transponder (PIT) tag reader (RePTaR) for use in obtaining individual identity scans of brown treesnakes (*Boiga irregularis*) on Gu&aring;han. We investigated if this setup could help remotely scan individuals in a laboratory setting with hopes of eventual use of RePTaR in the field to augment monitoring work being done in targeted study areas. Data came in 3 forms and from 3 sources: 1) scans (where snakes were successfully read by the reader) downloaded from RePTaR microSD cards, 2) scanning success (1 = scan, 0 = no scan) and trial information (e.g., approximated distance from antenna) from processed video recordings from each trial (while snakes could scan from farther away than 2 inches, we recorded failure to scan only within this 2-in area around all antennas), and 3) individual snake trait or tagging information from processing captured snakes. We used this information to run logistic regressions investigating what covariates could help describe variation in scanning success.

# Table of Contents
Data folder - Please download from ScienceBase, see text file for details. The data archive contains individual snake trait information (SnakeSizes.csv), behavioral trial data collected from videos (RePTaR_trials_Aug2021_AllData_2.csv), and RePTaR scans from microSD cards (202185_829_stat1.csv [side 1 of experimental arena] and 202185_829_stat2.csv [side 2 of experimental arena])

Figures folder - contains output from figure prepartory scripts (ScanCounts.png, ScanDist.png, and SnakeInfo.png)

Models folder - contains JAGS code to run all models

Results folder - contains output from JAGS models for all models

Scripts folder - all code to prepare, format, analyze (LogisticRegressionScanningSuccess.R), and plot (PlottingPatternsRePTaR.R)

# Required Packages and Versions Used
tidyr_1.1.3

dplyr_1.0.7

lubridate_1.7.10

tidyverse_1.3.1

ggplot2_3.3.5

LaCroixColoR_0.1.0

jagsUI_1.5.2

# Details of Article
Details of this work can be found in the published journal article on this topic:

Amburgey SA, Prakash A, Yackel Adams AA, Siers S, Converse SJ (2025). Development and evaluation of the remote passive integrated transponder tag reader for customizable monitoring of wildlife. Wildlife Society Bulletin e1569.

# How to Use This Repository
Start by downloading Data from ScienceBase (text file details where to find this archive). Go to file in Scripts to summarize, prepare, and plot data; Results and figures will save to other folders and can be viewed there
