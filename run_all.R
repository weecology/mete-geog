## To recreate the analysis of geographical patterns of S and N and their
## influence on rarity

## dependencies:
## sqlite3
## 

## download, install, and update scripts for the ecological data retriever
## http://ecodataretriever.org/
##

## Data Processing -------------------------------------------------------------
## download datasets
library(ecoretriever)

dir.create('./log_files')
public_data = c('MCDB', 'Gentry', 'BBS', 'FIA')

for(i in seq_along(public_data)) {
    ecoretriever::install(public_data[i], 'mysql', log_dir = './log_files') 
}

## query datasets


## prepare environmental GIS layers

## query GIS layers at site coordinates


## Data Analysis ---------------------------------------------------------------
## model constraints

## cross validate and model selection

## Summarize Results -----------------------------------------------------------