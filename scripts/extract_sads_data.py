"""Downloads and extracts the raw data for METE SADs analyses

Downloads BBS, FIA, Gentry, and MCDB (NABC and CBC are not public)
Queries databases to extract the raw data for analysis

WARNING: These databases are large and take a long time to download and a long
time to query.

"""

import os

import getpass
import shutil

def get_raw_data(queries, engine='sqlite', file_name='downloaded_data.sqlite',
                 host=None, port=3306, user='root'):
    """Function to connect to database, create query tables, and output to CSV.
    """
    assert engine in ('sqlite', 'mysql'), 'Engine must be either sqlite or mysql'
    if engine == 'sqlite':
        import sqlite3 as sqlite_dbapi
        connection = sqlite_dbapi.connect(file_name)
    elif engine == 'mysql':
        import MySQLdb as dbapi
        p=getpass.getpass('Enter MySQL Password')
        connection = dbapi.connect(host=host, port=port, user=user, passwd=p)
    
    cursor = connection.cursor()

    for query in queries:
        cursor.execute("""DROP DATABASE IF EXISTS queries;""")
        cursor.execute("""CREATE DATABASE queries;""")
        cursor.execute(query)

    connection.commit()

# BBS

bbs_queries = [# Step 1. TAXONOMY_BIRDS TAXONOMY.AOU_IN linked to BBS counts.Aou
               # Left Join = ALL records from BBS counts are included 
               # and only values from TAXONOMY that match
               """
               CREATE TABLE queries.aous (aou INT(11)) 
               SELECT counts.Aou FROM BBS.counts 
               GROUP BY counts.Aou;
               """,
               """
               CREATE TABLE queries.aou_too_1 (aou INT(11)) 
               SELECT Aou AS AOU, TAXON_ORDER_OUT AS TOO FROM queries.aous 
               LEFT JOIN TAXONOMY_BIRDS.TAXONOMY 
               ON aous.Aou = TAXONOMY.AOU_IN 
               WHERE TAXONOMY.DIURNALLANDBIRD = 1;
               """,
               # Step 2. Group By AOU IN and by TOO WHERE DIURNAL LANDBIRD = 1 - 
               # to yield AOU_TOO                        
               """
               CREATE TABLE queries.aou_too (aou INT(11)) 
               SELECT AOU, TOO FROM queries.aou_too_1 
               GROUP BY AOU, TOO;
               """,
               # Step 3. Create table with SiteID, Year, RunType=1 from weather table
               """
               CREATE TABLE queries.weather_subquery
               SELECT (weather.statenum * 1000 + weather.Route) AS SiteID,
               weather.Year, weather.RunType
               FROM BBS.weather
               WHERE weather.RunType = 1 AND weather.RPID = 101;
               """,
               # Step 4. To create table with SiteID, Year, Sp, abund:
               # Link together AOU_TOO and BBS Counts by AOU - 
               # Group By SiteID = state*1000 + route - TOO - Year 
               # Sum SpeciesTotal
               """
               CREATE TABLE queries.counts_too
               SELECT (counts.statenum * 1000) + counts.Route AS SiteID,
               counts.Year, aou_too.TOO, counts.RPID,
               SUM(counts.SpeciesTotal) AS AB 
               FROM BBS.counts INNER JOIN queries.aou_too ON 
               counts.Aou = aou_too.AOU
               GROUP BY SiteID, counts.Year, aou_too.TOO, counts.RPID
               HAVING (((counts.Year = 2009) AND (counts.RPID = 101)));
               """,
               """
               SELECT counts_too.SiteID, counts_too.Year, counts_too.TOO, 
               counts_too.AB
               FROM queries.counts_too INNER JOIN queries.weather_subquery
               ON counts_too.SiteID = weather_subquery.SiteID 
               AND counts_too.Year = weather_subquery.Year
               INTO OUTFILE '/tmp/bbs_spab.csv'
               FIELDS TERMINATED BY ',' 
               LINES TERMINATED BY '\n';
               """]

# CBC

cbc_queries = [# Step 1. Group By SPECIES_CODE and by TOO WHERE DIURNAL LANDBIRD = 1 - 
               # to yield sp_too
               """
               CREATE TABLE queries.spcodes
               SELECT SPECIES.SPECIES_CODE FROM CBC.SPECIES 
               GROUP BY SPECIES.SPECIES_CODE ;
               """, 
               """
               CREATE TABLE queries.sp_too_1
               SELECT SPECIES_CODE AS SPCODE, TAXON_ORDER_OUT AS TOO FROM queries.spcodes 
               LEFT JOIN TAXONOMY_BIRDS.TAXONOMY 
               ON spcodes.SPECIES_CODE = TAXONOMY.CBCSPCODE
               WHERE TAXONOMY.DIURNALLANDBIRD = 1 ;
               """,
               """
               CREATE TABLE queries.sp_too 
               SELECT SPCODE, TOO FROM queries.sp_too_1 
               GROUP BY SPCODE, TOO ;
               """,
               # Step 2. LINK OBS and SUB_AUX by SUB_ID, adding COUNT_YR = 109 (2008-2009)
               # to yield OBSDATA_CTYR_STEP1 (SUB_ID, COUNT_YR, SPECIES_CODE, HOW_MANY)               
               """
               CREATE TABLE queries.obs_1
               SELECT SUB_AUX.SUB_ID, SUB_AUX.COUNT_YR, OBS.SPECIES_CODE, OBS.HOW_MANY
               FROM CBC.SUB_AUX INNER JOIN CBC.OBS ON 
               SUB_AUX.SUB_ID = OBS.SUB_ID
               WHERE SUB_AUX.COUNT_YR = 109 ;
               """,
               # Step 3. LINK obs_1 to SUB by SUB_ID to add LOC_ID
               """
               CREATE TABLE queries.obs_2
               SELECT SUB.LOC_ID, obs_1.* 
               FROM queries.obs_1 INNER JOIN CBC.SUB ON 
               obs_1.SUB_ID = SUB.SUB_ID ;
               """,
               # Step 4. LINK obs_2 to LOC by LOC_ID to remove records WHERE SUBNATIONAL1_CODE 
               # IS NOT 'US-HI' OR 'US-AK', AND COUNTRY CODE IS 'CA', 'US', OR 'US-CA'
               """
               CREATE TABLE queries.obs_3
               SELECT obs_2.* 
               FROM queries.obs_2 INNER JOIN CBC.LOC ON
               obs_2.LOC_ID = LOC.LOC_ID
               WHERE LOC.COUNTRY_CODE NOT IN ('PR', 'DO', 'BS', 'BM', 'VI') AND 
               LOC.SUBNATIONAL1_CODE NOT IN ('US-HI') ;
               #WHERE LOC.COUNTRY_CODE IN ('CA', 'US', 'MX', 'US-CA') AND 
               #LOC.SUBNATIONAL1_CODE NOT IN ('US-HI', 'US-AK') ;
               """,
               # Step 5. CREATE obs_4 by adding sp_too.TOO
               """
               CREATE TABLE queries.obs_4
               SELECT obs_3.*, sp_too.TOO
               FROM queries.obs_3 INNER JOIN queries.sp_too 
               ON obs_3.SPECIES_CODE = sp_too.SPCODE ;
               """,
               # Step 6. To create table with Loc_id - Year - Duration_hrs - Num_Obs - Sp - abund:
               # GROUPing BY LOC_ID and SPECIES_CODE and SUMmming over HOW_MANY FROM obs_4
               """
               CREATE TABLE queries.obs_5
               SELECT obs_4.LOC_ID, obs_4.COUNT_YR AS YEAR, obs_4.TOO, 
               SUM(obs_4.HOW_MANY) AS AB
               FROM queries.obs_4 
               GROUP BY obs_4.LOC_ID, obs_4.COUNT_YR, obs_4.TOO ;
               """,
               # Step 7. Create table with the unique site ids
               """
               CREATE TABLE queries.loc_id
               SELECT DISTINCT LOC_ID FROM queries.obs_5 ; 
               """,
               # Step 8. Remove zeroes from data, and save table to file
               """
               SELECT 'site_id', 'year', 'too', 'ab'
               UNION ALL 
               SELECT * FROM queries.obs_5 
               WHERE obs_5.AB > 0
               INTO OUTFILE '/tmp/cbc_spab.csv'
               FIELDS TERMINATED BY ',' 
               LINES TERMINATED BY '\n' ;
               """,
               # Step 9. Create a table that has sampling effort measured as duration
               # for COUNT_YR = 109 (2008-2009)
               """
               CREATE TABLE queries.effort
               SELECT SUB.LOC_ID, SUB.SUB_ID, SUB.DURATION_HRS, SUB_AUX.COUNT_YR 
               FROM CBC.SUB INNER JOIN CBC.SUB_AUX
               ON SUB.SUB_ID = SUB_AUX.SUB_ID
               WHERE SUB_AUX.COUNT_YR = 109 ;
               """,
               # Step 10. Export coordinates and measures of sampling effort for only
               # those sites included in the species query
               """
               SELECT 'site_id', 'year', 'duration_hrs', 'latitude', 'longitude'
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
               """]

# FIA
fia_queries = [# Step 1: Select CN field from SURVEY table where ANN_INVENTORY = 'Y'
               """
               CREATE TABLE fia_queries.survey1
               SELECT cn FROM FIA.SURVEY
               WHERE SURVEY.ann_inventory = 'Y' ; 
               """, 
               # Step 2: Limit PLOT table to only some of the fields and only certain plots
               """
               CREATE TABLE fia_queries.plot1
               SELECT cn, srv_cn, cty_cn, invyr, statecd, unitcd,
                      countycd, plot, lat, lon, elev
               FROM FIA.PLOT
               WHERE 
                   # sampled at least one accessible forest land condition
                   PLOT.plot_status_cd = 1 AND    
                   # exclude periodic inventory plot (kindcd=0) and 
                   # eclude modeled periodic inventory (kindcd=4)
                   PLOT.kindcd NOT IN (0, 4) AND  
                   # all basically sampled the same way 
                   # (see Appendix B of Phase 2 documentation for explanations)
                   PLOT.designcd IN (1, 311, 312, 313, 314, 328, 220, 240) AND 
                   # limits data to those collected using the standardized methodology
                   # described in the National FIA Field Guide
                   PLOT.manual >= 1 AND
                   # Standard production plot 1 OR hot check 7 
                   # hot check just indicates that an expert was present at sampling
                   PLOT.qa_status IN (1, 7) AND 
                   # Field visited (as opposed to 2 = remote sensed)
                   PLOT.samp_method_cd = 1 AND
                   # drop Western Phase 3 plots that are'off subpanel'
                   PLOT.invyr != 9999 ; 
              """,
              # Step 3: Join these first two tables (fia_queries.fia_survey1 and fia_plot1)
              # to subset out only the plot records where ANN_INVENTORY is Yes
              """
              CREATE TABLE fia_queries.plot2
              SELECT plot1.* FROM fia_queries.plot1
              INNER JOIN fia_queries.survey1 ON
              plot1.srv_cn = survey1.cn 
              ORDER BY plot1.statecd, plot1.unitcd, plot1.countycd, plot1.plot, 
              plot1.invyr 
              DESC ;
              """,
              # Step 4: Use fia_queries.fia_plot2 to then generate intermediate table that
              # creates a unique id for each plot and includes the maximum value of
              # INVYR (inventory year)
              # the combination of STATECD, UNITCD, COUNTYCD, and PLOT can be used to 
              # generate a unique plot id
              """
              CREATE TABLE fia_queries.plot3
              SELECT cn, statecd, unitcd, countycd, plot, 
              MAX(invyr) AS invyr,
              AVG(lat) AS lat,
              AVG(lon) AS lon,
              AVG(elev) AS elev
              FROM fia_queries.plot2
              GROUP BY statecd, unitcd, countycd, plot ;
              """,
              # Step 5: Select PLT_CN field from COND table WHERE:
              """
              CREATE TABLE fia_queries.cond1
              SELECT plt_cn FROM FIA.COND 
              WHERE 
                  # Stand origin code (stdorgcd) is = 0: Natural stands
                  (stdorgcd = 0 OR stdorgcd is NULL) AND
                  # No observable treatment (i.e., no cutting, girdling, herbicides, etc.)
                  (trtcd1 = 0 OR trtcd1 is NULL) AND
                  (trtcd2 = 0 OR trtcd2 is NULL) AND
                  (trtcd3 = 0 OR trtcd3 is NULL) 
              GROUP BY plt_cn ;
              """,
              # Step 6: Use fia_cond1 to further subset fia_plot3
              # query took somewhere between 8 and 15 hours to run 
              # 136304 records
              """
              CREATE TABLE fia_queries.plot4
              SELECT plot3.* FROM fia_queries.plot3
              INNER JOIN fia_queries.cond1 ON plot3.cn = cond1.plt_cn ;
              """,
              # Step 7: subset the TREE table (14,817,974 records)
              """
              CREATE TABLE fia_queries.tree1
              SELECT TREE.PLT_CN, TREE.STATECD, TREE.UNITCD, TREE.COUNTYCD, 
                     TREE.PLOT, TREE.SPCD FROM FIA.TREE
              WHERE TREE.STATUSCD = 1 ;
              """,
              # Step 8: Use fia_plot4 to subset the fia_tree1 table using
              # CN = PLT_CN, STATECD, UNITCD, COUNTYCD, PLOT, INVYR
              """
              CREATE TABLE fia_queries.tree2
              SELECT ((tree1.statecd * 10000000000) +
                      (tree1.unitcd * 1000000000) +
                      (tree1.countycd * 1000000) +
                      tree1.plot) AS PlotID, 
                      tree1.plt_cn, tree1.spcd
              FROM fia_queries.tree1
              INNER JOIN fia_queries.plot4 ON
              tree1.plt_cn = plot4.cn 
              WHERE plot4.invyr < 3000 ; # change to != 9999
              """,
              # Step 9: Calculate totals per species per plot from Step 8 results
              """
              CREATE TABLE fia_queries.tree3
              SELECT tree2.PlotID, tree2.plt_cn, tree2.spcd, COUNT(tree2.spcd) AS AB
              FROM fia_queries.tree2 
              GROUP BY tree2.PlotID, tree2.PLT_CN, tree2.SPCD ;
              """,
              # Step 10: export species abu table
              """
              SELECT 'SiteID', 'plt_cn', 'spID', 'ab' 
              UNION ALL
              SELECT * FROM fia_queries.tree3
              INTO OUTFILE '/tmp/fia_spab.csv'
              FIELDS TERMINATED BY ',' 
              LINES TERMINATED BY '\n' ;
              """,
              # Step 11: create table of unique tree cn ids
              """
              CREATE TABLE fia_queries.tree_plots
              SELECT DISTINCT plt_cn AS plt_cn FROM fia_queries.tree3 ;
              """,
              # Step 12: export site coords table
              """
              # Takes approximately 10 minutes
              SELECT 'SiteID', 'cn', 'statecd', 'unitcd', 'countycd', 'plot',
                     'invyr', 'latitude', 'longitude', 'elev'
              UNION ALL
              SELECT ((statecd * 10000000000) +
                      (unitcd * 1000000000) +
                      (countycd * 1000000) +
                      plot) AS PlotID, plot4.*
              FROM fia_queries.plot4
              INNER JOIN fia_queries.tree_plots ON 
              tree_plots.plt_cn = plot4.cn
              INTO OUTFILE '/tmp/fia_coords.csv'
              FIELDS TERMINATED BY ',' 
              LINES TERMINATED BY '\n' ;
              """]

gentry_queries = ["""
                  # Export a species count table
                  SELECT 'siteID', 'spID', 'spName', 'ab' 
                  UNION ALL
                  SELECT counts.site_code, counts.species_id, 
                  CONCAT(species.genus, " ", species.species) AS species_name,
                  counts.count AS ab
                  FROM Gentry.counts INNER JOIN Gentry.species ON
                  counts.species_id = species.species_id
                  WHERE species.full_id = 1 
                  INTO OUTFILE '/tmp/gentry_spab.csv'
                  FIELDS TERMINATED BY ',' 
                  LINES TERMINATED BY '\n' ;
                  """, 
                  """                  
                  # Export coordinates file
                  SELECT 'siteID', 'latitude', 'longitude' 
                  UNION ALL
                  SELECT abbreviation, lat, lon
                  FROM Gentry.sites
                  INTO OUTFILE '/tmp/gentry_coords.csv'
                  FIELDS TERMINATED BY ',' 
                  LINES TERMINATED BY '\n' ;
                  """]

# NABC
nabc_queries = ["""CREATE TABLE queries.nabc_sample_dates 
                   SELECT NABA_2009.Count_State, NABA_2009.Count_Name, NABA_2009.Date
                   FROM NABC.NABA_2009
                   GROUP BY NABA_2009.Count_State, NABA_2009.Count_Name,
                   NABA_2009.Date;""",
                """CREATE TABLE queries.nabc_min_sample_date
                   SELECT nabc_sample_dates.Count_State,
                   nabc_sample_dates.Count_Name, 
                   MIN(nabc_sample_dates.Date) AS MinDate
                   FROM queries.nabc_sample_dates
                   GROUP BY nabc_sample_dates.Count_State, 
                   nabc_sample_dates.Count_Name;""",
                """CREATE TABLE queries.nabc_sp_ab_2009a
                   SELECT NABA_2009.Count_State, NABA_2009.Count_Name, 
                   CONCAT(NABA_2009.Count_State,"_",Left(NABA_2009.Count_Name,5),
                   "_",Right(NABA_2009.Count_Name,5)) AS SiteID,
                   YEAR(NABA_2009.Date) AS Year, NABA_2009.Scientific_Name, 
                   Sum(NABA_2009.Number_Butterflies) AS AB
                   FROM NABC.NABA_2009 INNER JOIN queries.nabc_min_sample_date ON 
                   (nabc_min_sample_date.Count_State = NABA_2009.Count_State) AND 
                   (nabc_min_sample_date.Count_Name = NABA_2009.Count_Name) AND 
                   (nabc_min_sample_date.MinDate = NABA_2009.Date)
                   GROUP BY NABA_2009.Count_State, NABA_2009.Count_Name,
                   NABA_2009.Date, NABA_2009.Scientific_Name;""",
                """CREATE TABLE queries.nabc_sp_ab_2009
                   SELECT nabc_sp_ab_2009a.SiteID, nabc_sp_ab_2009a.Year,
                   CONCAT(NABA_species.Genus, "_", NABA_species.Species) AS SpID, 
                   SUM(nabc_sp_ab_2009a.AB) AS Abund
                   FROM queries.nabc_sp_ab_2009a INNER JOIN NABC.NABA_species ON 
                   NABA_species.Scientific_Name = nabc_sp_ab_2009a.Scientific_Name 
                   GROUP BY nabc_sp_ab_2009a.SiteID, NABA_species.Genus, 
                   NABA_species.Species;""",
                # Dump into csv file, removing two sites that include significant 
                # outliers in abundance, one being the NABA Butterfly park
                """SELECT nabc_sp_ab_2009.* FROM queries.nabc_sp_ab_2009
                   WHERE nabc_sp_ab_2009.SiteID != "TX_NABA _ Park" AND 
                   nabc_sp_ab_2009.SiteID != "MN_Bear _ction"
                   INTO OUTFILE '/tmp/naba_spab.csv'
                   FIELDS TERMINATED BY ',' 
                   LINES TERMINATED BY '\n';"""]


#workdir = os.getcwd()
#shutil.copy('/tmp/bbs_spab.csv', '/home/kate/data/bbs_spab.csv')
#shutil.copy('/tmp/cbc_spab.csv', '/home/kate/data/cbc_spab.csv')
#shutil.copy('/tmp/fia_spab.csv', '/home/kate/data/fia_spab.csv')
#shutil.copy('/tmp/gentry_spab.csv', '/home/kate/data/gentry_spab.csv')
#shutil.copy('/tmp/mcdb_spab.csv', '/home/kate/data/mcdb_spab.csv')
#shutil.copy('/tmp/naba_spab.csv', '/home/kate/data/naba_spab.csv')