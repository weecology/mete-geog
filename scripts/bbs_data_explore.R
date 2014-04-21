## Explore the BBS dataset ---------------------------------------------

library(maps)

bbs08 = read.csv('./data/raw/bbs_2008_2012_spab.csv', header=F)
bbs12 = read.csv('./data/raw/bbs_2012_spab.csv', header=F)
coords = read.csv('./data/raw/bbs_coords.csv', header=F)

pdf('./figs/bbs_map.pdf')
map('world', c('canada', 'usa','mexico'),
    xlim=c(-135, -50), ylim=c(25, 60))
map('state', add=T)
true = coords[,1] %in% bbs12[,1]
points(coords[true, 3], coords[true, 2], cex=.5)
true = coords[,1] %in% bbs08[,1]
points(coords[true, 3], coords[true, 2], cex=.25,
       pch=19, col='red')
legend('bottomright', c('2012', '2008-12'),
       pch=c(1, 19), col=c('black','red'), bty='n')
dev.off()

## BBS MySQL queries ----------------------------------------------------------

# Step 1. TAXONOMY_BIRDS TAXONOMY.AOU_IN linked to BBS counts.Aou
# Left Join = ALL records from BBS counts are included 
# and only values from TAXONOMY that match
# Step 2. Group By AOU IN and by TOO WHERE DIURNAL LANDBIRD = 1 - 
# to yield AOU_TOO

DROP DATABASE IF EXISTS queries;
CREATE DATABASE queries;

CREATE TABLE queries.aous (aou INT(11)) 
SELECT counts.Aou FROM BBS.counts 
GROUP BY counts.Aou ;

CREATE TABLE queries.aou_too_1 (aou INT(11)) 
SELECT Aou AS AOU, TAXON_ORDER_OUT AS TOO FROM queries.aous 
LEFT JOIN TAXONOMY_BIRDS.TAXONOMY 
ON aous.Aou = TAXONOMY.AOU_IN 
WHERE TAXONOMY.DIURNALLANDBIRD = 1 ;

CREATE TABLE queries.aou_too (aou INT(11)) 
SELECT AOU, TOO FROM queries.aou_too_1 
GROUP BY AOU, TOO ;

# 3. Create table with SiteID, Year, RunType=1 from weather table
CREATE TABLE queries.weather_sub
SELECT (weather.statenum * 1000 + weather.Route) AS SiteID,
weather.Year, weather.RunType
FROM BBS.weather
WHERE weather.RunType = 1 AND weather.RPID = 101 AND weather.Year >= 2008;

# 4. Create table to record how many times from 2008 to 2012 (5 years)
# that a particular site was sampled. Only sites that were sampled
# continuously for the past 5 years will be used in the 5 yr window
# dataset.

CREATE TABLE queries.weather_cnt
SELECT SiteID, COUNT(SiteID) AS Count 
FROM queries.weather_sub
GROUP BY SiteID;

# 5. distribute the year count information into the year specific 
# weather records 

CREATE TABLE queries.weather_sub_cnt
SELECT weather_sub.SiteID, weather_sub.Year, weather_sub.RunType,
weather_cnt.Count
FROM queries.weather_sub INNER JOIN queries.weather_cnt
ON weather_sub.SiteID = weather_cnt.SiteID ;


# 5. To create table with SiteID, Year, Sp, abund:
# Link together AOU_TOO and BBS Counts by AOU - 
# Group By SiteID = state*1000 + route - TOO - Year 
# Sum SpeciesTotal for 2012

CREATE TABLE queries.counts_too
SELECT (counts.statenum * 1000) + counts.Route AS SiteID,
counts.Year, aou_too.TOO, counts.RPID,
SUM(counts.SpeciesTotal) AS AB 
FROM BBS.counts INNER JOIN queries.aou_too
ON counts.Aou = aou_too.AOU
GROUP BY SiteID, counts.Year, aou_too.TOO, counts.RPID
HAVING (((counts.Year >= 2008) AND (counts.RPID = 101))) ;

# 6. Output results for 2012 contiguous sampling

SELECT counts_too.SiteID, counts_too.Year, counts_too.TOO, 
counts_too.AB
FROM queries.counts_too INNER JOIN queries.weather_sub_cnt
ON counts_too.SiteID = weather_sub_cnt.SiteID 
AND counts_too.Year = weather_sub_cnt.Year
AND counts_too.Year = 2012
INTO OUTFILE '/tmp/bbs_2012_spab.csv'
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n';  

# 7. Output results for the 2008 to 2012 contiguous sampling

SELECT counts_too.SiteID, counts_too.TOO, SUM(counts_too.AB) AS AB
FROM queries.counts_too INNER JOIN queries.weather_sub_cnt
ON counts_too.SiteID = weather_sub_cnt.SiteID
AND counts_too.Year = weather_sub_cnt.Year
AND weather_sub_cnt.Count = 5
GROUP BY counts_too.SiteID, counts_too.TOO
INTO OUTFILE '/tmp/bbs_2008_2012_spab.csv'
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n';  

# 8. Output coordinates for routes examined:

SELECT (routes.statenum * 1000) + routes.route AS SiteID, lati, loni
FROM BBS.routes 
INTO OUTFILE '/tmp/bbs_coords.csv'
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n';  
