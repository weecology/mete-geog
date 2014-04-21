## examine MCDB datastructure

library(ecoretriever)
library(maps)

MCDB = fetch('MCDB')

names(MCDB)

attach(MCDB)

## Errors detected in MCDB ---------------------------------------------

## trapping
# 3 records have a final_year of 200 instead of 2000
trapping_errors = subset(trapping, final_year == 200)
write.csv(trapping_errors, file ='./data/mcdb_trapping_errors.csv')


## why doesn't the variable n_sampling_months always correspond with the 
## variable months_of_sampling ?
## for example:
site_id_ex = sites$site_id[sites$time_series==1][5]
trapping[trapping$site_id == site_id_ex, ]

## sites
# Longitude is misspelled
# two sites have longitudes are too negative
summary(sites$lonitude)
sum(sites$lonitude < -180, na.rm=T)

## some of the study durations are shorter than the number of years
sum(sites$study_duration < sites$n_years, na.rm=T)

sites_errors = subset(sites, lonitude <= -180 | study_duration < n_years)
write.csv(sites_errors, file='./data/mcdb_sites_errors.csv')

## Fix errors in MCDB -----------------------------
trapping$final_year[trapping$final_year == 200] = 2000
sites$n_years[sites$study_duration < sites$n_years] = 1
names(sites) = c(names(sites)[1:7], 'longitude', names(sites)[9:20])
sites$longitude = ifelse(sites$longitude < -180, 
                         sites$longitude * 1e-6,
                         sites$longitude)
                        

## Explore data ----------------------------------------


years = unlist(sapply(1:nrow(trapping), function(i)
               if (!is.na(trapping$initial_year[i]) & !is.na(trapping$final_year[i]))
                 trapping$initial_year[i]:trapping$final_year[i]))

site_ids = unlist(sapply(1:nrow(trapping), function(i)
                  if (!is.na(trapping$initial_year[i]) & !is.na(trapping$final_year[i]))
                    rep(trapping$site_id[i],
                        length(trapping$initial_year[i]:trapping$final_year[i]))))


pdf('./figs/mcdb_years.pdf', width=7*2)
par(mfrow=c(1,2))
hist(years)
abline(v=median(years), col='red')
hist(years[years > 1980])
abline(v=median(years[years>1980]), col='red')
dev.off()

median_year = median(years[years>1980])
median_year
#[1] 1997

## median is 1997 given a 10 year time span
summary((median_year - 5) : (median_year + 5))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1992    1994    1997    1997    2000    2002
gd_years = years >= 1992 & years <= 2002
length(years[gd_years])
# 611 sampling events during that period of time
uni_site_ids = unique(site_ids[gd_years])
length(uni_site_ids)
# 353 unique sites sampled 
tab = table(sites$time_series[sites$site_id %in% uni_site_ids])
tab
#   0   1 
# 288  65 
tab = table(sites$study_duration[sites$site_id %in% uni_site_ids])
n_months = as.numeric(names(tab))
plot(n_months, cumsum(tab))
cumsum(tab[n_months <= 12])
#  1   2   3   4   5   7   8   9  10  11  12 
# 49  93 137 160 168 171 172 173 177 178 192 

## define good site ids with following criteria:
## 1. study duration is less than 1 year
## 2. coordinates exist
## 3. abundance records are 'raw' not 'estimated'

true  = sites$site_id %in% uni_site_ids &
        sites$study_duration <= 12 &
        !is.na(sites$longitude) &
        !is.na(sites$latitude) &
        sites$abundance_data_format == 'raw'
        
gd_sites = sites$site_id[true]

## now use these site ids to pull out good sampling events from trapping
## only want a single sampling event per site
## so first pull all data then subset down to the event closest to median year
trapping = trapping[trapping$site_id %in% gd_sites, ]

trap_uni = as.data.frame(matrix(NA, nrow=0, ncol=ncol(trapping)))
names(trap_uni) = names(trapping)
for (i in seq_along(gd_sites)) {
    tmp_trap = trapping[trapping$site_id == gd_sites[i], ]
    closest_yr = min(tmp_trap$final_year - median_year)
    index = which(tmp_trap$final_year - median_year == closest_yr)
    trap_uni  = rbind(trap_uni , tmp_trap[index, ])
}

nrow(trap_uni)
#[1] 192

## define good trapping events with following criteria:
## 1. has trap_nights information recorded
## 2. does not include pitfall traps
## 3. does not include large traps
## 4. can include any combination or small and snap traps
true = !is.na(trap_uni$trap_nights) &
       trap_uni$pitfall_traps == 0 &
       trap_uni$large_traps == 0 

gd_traps = trap_uni[true, ]
nrow(gd_traps)
#[1] 115

gd_trap_sites = gd_traps$site_id

pdf('./figs/mcdb_map.pdf')
map('world')
mtext(side=3, 'All sites')
points(sites$longitude, sites$latitude, pch=1, cex=.5, col='red')
true = sites$site_id %in% gd_trap_sites
points(sites$longitude[true], sites$latitude[true], pch=19, cex=.25, col='blue')
##
map('world')
mtext(side=3, 'Good sites')
points(sites$longitude[true], sites$latitude[true], pch=1, cex=.5, col='red')
dev.off()

boxplot(sites$n_years)

table(sites$time_series)
#  0   1 
#921  79
## so 921 sites we have data aggregated across sampling events

## how many of these sites with aggregation have multiple sampling events
tab = table(sites$study_duration[sites$time_series==0])
plot(as.numeric(names(tab)), cumsum(tab), log='x', type='o')
abline(v=12)
##
tab = table(sites$study_duration[sites$time_series==1])
plot(as.numeric(names(tab)), cumsum(tab), log='x', type='o')
abline(v=12)

## hypothetical example
site_id_ex = sites$site_id[sites$time_series==1][2]
trapping[trapping$site_id == site_id_ex, ]


par(mfrow=c(1,2))
hist(sites$n_years[sites$time_series == 0])
hist(sites$n_years[sites$time_series == 1])

summary(sites$n_years[sites$time_series == 0])

plot(sites$study_duration, sites$n_years)
abline(a=0,b=1)

par(mfrow=c(1,2))
plot(trapping$n_sampling_months, log10(trapping$trap_nights))
plot(trapping$n_sampling_months, log10(trapping$trap_nights), xlim=c(0, 12))

plot(trapping$initial_year

hist(trapping$initial_year)
hist(trapping$final_year)

## MCDB MySQL Queries ---------------------------------------------------------
SHOW TABLES FROM queries ;

DROP DATABASE IF EXISTS queries;
CREATE DATABASE queries;

SELECT * FROM MCDB_trapping
WHERE small_traps = 1 AND
trap_nights is not null AND
pitfall_traps = 0 AND
large_traps = 0 AND
snap_traps = 0


mcdb_queries = ["""CREATE TABLE queries.mcdb_abu
                SELECT communities.Site_ID AS site,
                communities.Initial_year AS year,
                communities.Species_ID AS sp,
                communities.Abundance AS abu,
                sites.Abundance_data_format AS format 
                FROM MCDB.communities
                INNER JOIN MCDB.sites USING (Site_ID)
                HAVING (((abu > 0) AND (format = 'raw')));""",
                # Create intermediate table that has the value of the
                # earliest year of 
                # sampling for each site 
                """CREATE TABLE queries.mcdb_minyear
                SELECT mcdb_abu.site, Min(mcdb_abu.year) AS year
                FROM queries.mcdb_abu
                GROUP BY mcdb_abu.site;""",
                # Create intermidate table that 
                # Use intermediate table to select out only one year of data
                # per site
                """CREATE TABLE queries.mcdb_out
                SELECT mcdb_abu.site, mcdb_abu.year, mcdb_abu.sp, mcdb_abu.abu 
                FROM queries.mcdb_abu
                INNER JOIN queries.mcdb_minyear USING (site, year);""",
                # Dump into csv file
                """SELECT mcdb_out.* FROM queries.mcdb_out
                INTO OUTFILE '/tmp/mcdb_spab.csv'
                FIELDS TERMINATED BY ',' 
                LINES TERMINATED BY '\n';"""]


### miscellanous mysql commands
SHOW DATABASES;
USE MCDB ; 
SHOW TABLES ;
SELECT COLUMNS FROM communities ;
SELECT * FROM communities ;

SHOW TABLES FROM MCDB ;
SHOW TABLE STATUS FROM MCDB ;

SHOW COLUMNS FROM MCDB.communities ;

SHOW TABLES FROM queries ;
DROP DATABASE IF EXISTS queries;
SHOW TABLES FROM queries ;

CREATE TABLE queries.mcdb_abu
SELECT communities.Site_ID AS site,
communities.Initial_year AS year,
communities.Species_ID AS sp,
communities.Abundance AS ab,
sites.Abundance_data_format AS format
FROM MCDB.communities
INNER JOIN MCDB.sites USING (Site_ID)
HAVING (((ab > 0) AND (format = 'raw')));


SHOW TABLES FROM queries ;
SHOW COLUMNS FROM queries ;
SHOW COLUMNS FROM queries.mcdb_abu ;
SELECT DISTINCT sp FROM queries.mcdb_abu ;
SELECT * FROM queries.mcdb_abu ;
SELECT * FROM queries.mcdb_abu LIMIT 5;
SHOW COLUMNS FROM MCDB.communities ;

## SQL doesn't care about variable type setting
SELECT Site_ID FROM MCDB.communities ;
SELECT site_id FROM MCDB.communities ;
SELECT site_id FROM MCDB.communities ;

CREATE TABLE queries.mcdb_abu
SELECT communities.Site_ID AS site,
communities.Initial_year AS year,
communities.Species_ID AS sp,
communities.Abundance AS abu,
sites.Abundance_data_format AS format 
FROM MCDB.communities
INNER JOIN MCDB.sites USING (Site_ID)
HAVING (((abu > 0) AND (format = 'raw')));
CREATE TABLE queries.mcdb_abu                    SELECT communities.Site_ID AS site,                    communities.Initial_year AS year,                    communities.Species_ID AS sp,                    communities.Abundance AS abu,                    sites.Abundance_data_format AS format                     FROM MCDB.communities                    INNER JOIN MCDB.sites USING (Site_ID)                    HAVING (((abu > 0) AND (format = 'raw')));
CREATE TABLE queries.mcdb_minyear
SELECT mcdb_abu.site, Min(mcdb_abu.year) AS year
FROM queries.mcdb_abu
GROUP BY mcdb_abu.site;
CREATE TABLE queries.mcdb_minyear                    SELECT mcdb_abu.site, Min(mcdb_abu.year) AS year                    FROM queries.mcdb_abu                    GROUP BY mcdb_abu.site;
CREATE TABLE queries.mcdb_gdsites 
SELECT Site_ID FROM MCDB.trapping
WHERE trap_nights is not null AND
small_traps = 1 AND
pitfall_traps = 0 AND
large_traps = 0 AND
snap_traps = 0;
CREATE TABLE queries.mcdb_gdsites                     SELECT Site_ID FROM MCDB.trapping                    WHERE trap_nights is not null AND                          small_traps = 1 AND                          pitfall_traps = 0 AND                          large_traps = 0 AND                          snap_traps = 0;
CREATE TABLE queries.mcdb_gdsites 
SELECT Site_ID AS site 
FROM MCDB.trapping
WHERE trap_nights is not null AND
small_traps = 1 AND
pitfall_traps = 0 AND
large_traps = 0 AND
snap_traps = 0;
CREATE TABLE queries.mcdb_gdsites                     SELECT Site_ID AS site                     FROM MCDB.trapping                    WHERE trap_nights is not null AND                          small_traps = 1 AND                          pitfall_traps = 0 AND                          large_traps = 0 AND                          snap_traps = 0;
DROP TABLE queries.mcdb_gdsites
;
CREATE TABLE queries.mcdb_gdsites                     SELECT Site_ID AS site                     FROM MCDB.trapping                    WHERE trap_nights is not null AND                          small_traps = 1 AND                          pitfall_traps = 0 AND                          large_traps = 0 AND                          snap_traps = 0;
SELECT * FROM queries.mcdb_minyear
;
show tables queries.mcdb_gdsites
;
SHOW TABLES FROM queries ; 
SELECT * queries.mcdb_gdsites
;
SELECT * queries.mcdb_gdsites
;
SELECT * queries.mcdb_gdsites;
SELECT * queries.mcdb_minyear;
SELECT * FROM queries.mcdb_minyear;
SELECT DISTINCT site FROM queries.mcdb_minyear;
SELECT * FROM queries.mcdb_minyear;


