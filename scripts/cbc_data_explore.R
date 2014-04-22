library(maps)

spab = read.csv('./data/raw/cbc_spab.csv')
coords = read.csv('./data/raw/cbc_coords.csv')

pdf('./figs/cbc_map.pdf')
map('world', c('canada', 'usa','mexico'),
    xlim=c(-180, -50), ylim=c(25, 70))
map('state', add=T)
points(coords$longitude, coords$latitude, cex=.25,
       pch=19, col='red')
dev.off()

# CBC MySQL queries -----------------------------------------------------------
# Step 1. Group By SPECIES_CODE and by TOO WHERE DIURNAL LANDBIRD = 1 - 
# to yield sp_too
DROP DATABASE IF EXISTS queries;
CREATE DATABASE queries ;

CREATE TABLE queries.spcodes
SELECT SPECIES.SPECIES_CODE FROM CBC.SPECIES 
GROUP BY SPECIES.SPECIES_CODE ;

CREATE TABLE queries.sp_too_1
SELECT SPECIES_CODE AS SPCODE, 
TAXON_ORDER_OUT AS TOO FROM queries.spcodes 
LEFT JOIN TAXONOMY_BIRDS.TAXONOMY 
ON spcodes.SPECIES_CODE = TAXONOMY.CBCSPCODE
WHERE TAXONOMY.DIURNALLANDBIRD = 1 ;

CREATE TABLE queries.sp_too 
SELECT SPCODE, TOO FROM queries.sp_too_1 
GROUP BY SPCODE, TOO ;

# Step 2. LINK OBS and SUB_AUX by SUB_ID, adding COUNT_YR = 109 (2008-2009)
# to yield OBSDATA_CTYR_STEP1 (SUB_ID, COUNT_YR, SPECIES_CODE, HOW_MANY)

CREATE TABLE queries.obs_1
SELECT SUB_AUX.SUB_ID, SUB_AUX.COUNT_YR, OBS.SPECIES_CODE, OBS.HOW_MANY
FROM CBC.SUB_AUX INNER JOIN CBC.OBS ON 
SUB_AUX.SUB_ID = OBS.SUB_ID
WHERE SUB_AUX.COUNT_YR = 109 ;

# Step 3. LINK obs_1 to SUB by SUB_ID to add LOC_ID

CREATE TABLE queries.obs_2
SELECT SUB.LOC_ID, obs_1.* 
FROM queries.obs_1 INNER JOIN CBC.SUB ON 
obs_1.SUB_ID = SUB.SUB_ID ;

# Step 4. LINK obs_2 to LOC by LOC_ID to remove records WHERE SUBNATIONAL1_CODE 
# IS NOT 'US-HI' OR 'US-AK', AND COUNTRY CODE IS 'CA', 'US', OR 'US-CA'

CREATE TABLE queries.obs_3
SELECT obs_2.* 
FROM queries.obs_2 INNER JOIN CBC.LOC ON
obs_2.LOC_ID = LOC.LOC_ID
WHERE LOC.COUNTRY_CODE IN ('CA', 'US', 'US-CA') AND 
LOC.SUBNATIONAL1_CODE NOT IN ('US-HI', 'US-AK') ;

# Step 5. CREATE obs_4 by adding sp_too.TOO

CREATE TABLE queries.obs_4
SELECT obs_3.*, sp_too.TOO
FROM queries.obs_3 INNER JOIN queries.sp_too 
ON obs_3.SPECIES_CODE = sp_too.SPCODE ;

# Step 6. To create table with SiteID - Year - Duration_hrs - Num_Obs - Sp - abund:
# GROUPing BY LOC_ID and SPECIES_CODE and SUMmming over HOW_MANY FROM obs_4

CREATE TABLE queries.obs_5                
SELECT obs_4.LOC_ID, obs_4.COUNT_YR AS YEAR, obs_4.TOO, SUM(obs_4.HOW_MANY) AS AB
FROM queries.obs_4 
GROUP BY obs_4.LOC_ID, obs_4.COUNT_YR, obs_4.TOO ;

# Step 7. Create table with the unique site ids

CREATE TABLE queries.loc_id
SELECT DISTINCT LOC_ID FROM queries.obs_5 ; 

# Step 8. Remove zeroes from data, and save table to file

SELECT 'loc_id', 'year', 'too', 'ab'
UNION ALL 
SELECT * FROM queries.obs_5 
WHERE obs_5.AB > 0
INTO OUTFILE '/tmp/cbc_spab.csv'
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n' ;

# Step 9. Create a table that has sampling effort measured as duration
# for COUNT_YR = 109 (2008-2009)

CREATE TABLE queries.effort
SELECT SUB.LOC_ID, SUB.SUB_ID, SUB.DURATION_HRS, SUB_AUX.COUNT_YR 
FROM CBC.SUB INNER JOIN CBC.SUB_AUX
ON SUB.SUB_ID = SUB_AUX.SUB_ID
WHERE SUB_AUX.COUNT_YR = 109 ;

# Step 10. Export coordinates and measures of sampling effort for only those sites
# include in the species query

SELECT 'loc_id', 'duration_hrs', 'latitude', 'longitude'
UNION ALL 
SELECT effort.LOC_ID, effort.DURATION_HRS, LOC.LATITUDE, LOC.LONGITUDE
FROM queries.effort
INNER JOIN CBC.LOC
      ON effort.LOC_ID = LOC.LOC_ID
INNER JOIN queries.loc_id
      ON effort.LOC_ID = loc_id.LOC_ID
INTO OUTFILE '/tmp/cbc_coords.csv'
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n' ;


## examine and export  datafiles ----------------------------------------------
SHOW TABLES FROM CBC ;
show columns from CBC.AUX_LOC ;
show columns from CBC.AUX_WEATHER ;
show columns from CBC.COUNT_YEARS ;
show columns from CBC.LOC ;
show columns from CBC.OBS ;
show columns from CBC.SUB ;
show columns from CBC.SUB_AUX ;

SELECT * FROM CBC.SUB_AUX 
INTO OUTFILE '/tmp/cbc_sub_aux.csv'
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n';

SELECT * FROM CBC.SUB
INTO OUTFILE '/tmp/cbc_sub.csv'
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n';  

