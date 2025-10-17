import pymysql
import pylab as pl
import numpy as np
import pandas
import warnings
import scipy.stats as stats

import db_info
import Chi2Mfit_plotting

#ssh -N -L 12345:34.73.248.9:3306 iced@stoneage.ice-d.org

[myhost,myport,myuser,mypassword,mydatabase,mydatabase2] = db_info.credentials()

pandas.set_option("display.max_colwidth", None)
pandas.set_option('display.max_rows', None)
pandas.set_option('display.max_columns', None)
warnings.filterwarnings('ignore')

NHall = "(base_sample.lat_DD > 30)"
NPac = "(base_sample.lat_DD > 30 AND base_sample.lon_DD < -94)"
NAtl = "(base_sample.lat_DD > 30 AND base_sample.lon_DD > -94 AND base_sample.lon_DD < 59)"
HMA = "(base_sample.lat_DD > 30 AND base_sample.lon_DD > 59)"
no_HMA = "(base_sample.lat_DD > 30 AND base_sample.lon_DD < 59)"
region = NHall
tmin = 10000
tmax = 20000
age_scale = "LSDn"
site_count = 4
norm_camels = False
mc_runs = 1000

def reader_connect_to_db():
    dbc = pymysql.connect(host=myhost,port=myport,user=myuser,password=mypassword,database=mydatabase)

    return dbc

def bps_query():
    dbc = reader_connect_to_db()
    dbcursor = dbc.cursor()

    query = f"""SELECT base_sample.name, base_site.short_name, base_calculatedage.t_{age_scale}, base_calculatedage.dtint_{age_scale}, base_sample.lat_DD, base_sample.lon_DD, base_calculatedage.nuclide
        FROM base_sample
        JOIN base_site ON base_sample.site_id = base_site.id
        JOIN base_calculatedage ON base_calculatedage.sample_id = base_sample.id
        JOIN base_application_sites ON base_application_sites.site_id = base_site.id
        WHERE base_site.what LIKE "%moraine%"
        AND base_sample.site_id IN
        (SELECT base_sample.site_id FROM base_sample GROUP BY base_sample.site_id HAVING COUNT(*) >= {site_count})
        AND ({region} OR base_sample.lat_DD < -30)
        AND base_application_sites.application_id = 2
        AND base_calculatedage.nuclide LIKE "%N10quartz%"
        ORDER BY base_sample.lat_DD"""
    dbcursor.execute(query)
    result = dbcursor.fetchall()

    dbc.close()

    result_list = np.array(result)

    return result_list

def bps_queryALL():
    dbc = reader_connect_to_db()
    dbcursor = dbc.cursor()

    query = f"""SELECT base_sample.name, base_site.short_name, base_calculatedage.t_{age_scale}, base_calculatedage.dtint_{age_scale}, base_sample.lat_DD, base_sample.lon_DD, base_calculatedage.nuclide
        FROM base_sample
        JOIN base_site ON base_sample.site_id = base_site.id
        JOIN base_calculatedage ON base_calculatedage.sample_id = base_sample.id
        JOIN base_application_sites ON base_application_sites.site_id = base_site.id
        WHERE base_site.what LIKE "%moraine%"
        AND ({region} OR base_sample.lat_DD < -30)
        AND base_application_sites.application_id = 2
        AND (base_calculatedage.t_{age_scale} > {tmin - 1000} AND base_calculatedage.t_{age_scale} < {tmax + 1000})
        ORDER BY base_sample.lat_DD"""
    dbcursor.execute(query)
    result = dbcursor.fetchall()
    
    dbc.close()

    result_list = np.array(result)

    return result_list




def grouping(result_list):
    group_column = result_list[:,1]
    unique_groups = np.unique(group_column)
    
    grouped_data = {}

    for group_id in unique_groups:
        group_data = result_list[group_column == group_id]
        grouped_data[group_id] = group_data

    
    return grouped_data



def Chi2_filtering(group_result):

    master_array = []

    for group_id, group_data in group_result.items():
        ages = group_data[:,2].astype('float')
        errors = group_data[:,3].astype('float')
        errors[errors == 0] = 500
        site = group_id
        lat = group_data[:,4].astype('float')
        lon = group_data[:,5].astype('float')
        nuclide = group_data[:,6].astype('str')

        ages_init = ages
        errors_init = errors

        ages_final = ages
        errors_final = errors

        while True and len(ages_final) >= (len(ages_init) // 2):

            prune, ages_final, errors_final = Chi2Mfit_plotting.X2_pruning(ages_final, errors_final)

            if not prune:
                break
        
        if Chi2Mfit_plotting.norm_expected(ages_final, errors_final) > (tmin - 1000) and Chi2Mfit_plotting.norm_expected(ages_final, errors_final) < (tmax + 1000) and len(ages_final) > 2 and len(ages_final) >= (len(ages_init) // 2):
            
            site_array = np.array(list(zip(ages_final, errors_final, [site] * len(group_data), lat, lon, nuclide)), dtype=[('ages', 'float'), ('errors', 'float'), ('site', 'U10'), ('lat', 'float'), ('lon', 'float'), ('nuclide', 'U10')])
            master_array.append(site_array)


    return master_array



def data_organize(filtered_result):

    dataframes = []

    for array in filtered_result:
        df = pandas.DataFrame(array)
        dataframes.append(df)

    df1 = pandas.concat(dataframes, ignore_index=True)

    CriteriaNH = df1['lat'] > 0
    CriteriaSH = df1['lat'] < 0

    dfNH = df1[CriteriaNH].reset_index(drop=True)
    dfSH = df1[CriteriaSH].reset_index(drop=True)

    with open('/home/jtulenko/sw/analyses/bpseesaw_AGU23/dfNHfilter.txt', "w") as file1:
        file1.write(str(dfNH) + '\n')

    with open('/home/jtulenko/sw/analyses/bpseesaw_AGU23/dfSHfilter.txt', "w") as file2:
        file2.write(str(dfSH) + '\n')

    print('-----------------------')
    print(f'NH dataset filtered n = {len(dfNH)}')
    print(f'SH dataset filtered n = {len(dfSH)}')

    agesNH = dfNH['ages']
    errorsNH = dfNH['errors']

    agesSH = dfSH['ages']
    errorsSH = dfSH['errors']

    sum_agesNH = 0
    sum_agesSH = 0

    xi = np.arange(tmin,tmax,1)

    for e in range(0, len(agesNH)):
        ynh = (1/(errorsNH[e]*np.sqrt(2*np.pi)))*(np.exp(-0.5*(((xi-agesNH[e])/errorsNH[e])**2)))
        if norm_camels:
            ynh /= np.max(ynh)

        sum_agesNH = sum_agesNH + ynh

    for o in range(0, len(agesSH)):
        ysh = (1/(errorsSH[o]*np.sqrt(2*np.pi)))*(np.exp(-0.5*(((xi-agesSH[o])/errorsSH[o])**2)))
        if norm_camels:
            ysh /= np.max(ysh)

        sum_agesSH = sum_agesSH + ysh

    vlines = np.arange(tmin,tmax,500)

    
    pl.ion()
    pl.figure()
    pl.xlim(tmin,tmax)
    pl.ylim(0,(1.1 * np.max((sum_agesNH, sum_agesSH))))
    pl.xlabel('age [a]')
    pl.ylabel('Relative Probability')
    pl.yticks([])
    pl.title('Pseudo time-series summed pdfs FILTERED moraine ages NH (blue), SH (red)')
    pl.vlines(vlines, ymin=0, ymax=(1.1 * np.max((sum_agesNH, sum_agesSH))),linestyles='dashed')
    pl.plot(xi,sum_agesNH, 'b')
    pl.plot(xi, sum_agesSH, 'r')

    return sum_agesNH, sum_agesSH


def data_organizeALL(query_resultALL):
    df1 = pandas.DataFrame(list(query_resultALL))

    lat = df1[4].astype('float')

    CriteriaNH = lat > 0
    CriteriaSH = lat < 0

    dfNH = df1[CriteriaNH].reset_index(drop=True)
    dfSH = df1[CriteriaSH].reset_index(drop=True)

    with open('/home/jtulenko/sw/analyses/bpseesaw_AGU23/dfNHall.txt', "w") as file3:
        file3.write(str(dfNH) + '\n')

    with open('/home/jtulenko/sw/analyses/bpseesaw_AGU23/dfSHall.txt', "w") as file4:
        file4.write(str(dfSH) + '\n')

    print('-----------------------')
    print(f'NH dataset ALL n = {len(dfNH)}')
    print(f'SH dataset ALL n = {len(dfSH)}')

    agesNH = dfNH[2]
    errorsNH = dfNH[3]
    errorsNH[errorsNH == 0] = 500

    agesSH = dfSH[2]
    errorsSH = dfSH[3]
    errorsSH[errorsSH == 0] = 500

    sum_agesNH = 0
    sum_agesSH = 0

    xi = np.arange(tmin,tmax,1)

    for e in range(0, len(agesNH)):
        ynh = (1/(errorsNH[e]*np.sqrt(2*np.pi)))*(np.exp(-0.5*(((xi-agesNH[e])/errorsNH[e])**2)))
        if norm_camels:
            ynh /= np.max(ynh)

        sum_agesNH = sum_agesNH + ynh

    for o in range(0, len(agesSH)):
        ysh = (1/(errorsSH[o]*np.sqrt(2*np.pi)))*(np.exp(-0.5*(((xi-agesSH[o])/errorsSH[o])**2)))
        if norm_camels:
            ysh /= np.max(ysh)

        sum_agesSH = sum_agesSH + ysh

    vlines = np.arange(tmin,tmax,500)
    pl.ion()
    pl.figure()
    pl.xlim(tmin,tmax)
    pl.ylim(0,(1.1 * np.max((sum_agesNH, sum_agesSH))))
    pl.vlines(vlines, ymin=0, ymax=(1.1 * np.max((sum_agesNH, sum_agesSH))),linestyles='dashed')
    pl.xlabel('age [a]')
    pl.ylabel('Relative Probability')
    pl.yticks([])
    pl.title('Pseudo time-series summed pdfs ALL moraine ages NH (blue), SH (red)')
    pl.plot(xi,sum_agesNH, 'b')
    pl.plot(xi, sum_agesSH, 'r')

    return sum_agesNH, sum_agesSH

def data_organizeRAND(query_resultALL):
    df1 = pandas.DataFrame(list(query_resultALL))

    rand = np.random.randint(0, 10000, size=mc_runs)
    sum_agesu_list = []
    sum_agesi_list = []

    xi = np.arange(tmin,tmax,1)

    for i in range(0, len(rand)):

        shuffled_df = df1.sample(frac=1, random_state=rand[i])

        split_ratio = 0.75
        split_pt = int(len(shuffled_df) * split_ratio)

        dfu = shuffled_df.iloc[:split_pt].reset_index(drop=True)
        dfi = shuffled_df.iloc[split_pt:].reset_index(drop=True)

        agesu = dfu[2]
        errorsu = dfu[3]
        errorsu[errorsu == 0] = 500

        agesi = dfi[2]
        errorsi = dfi[3]
        errorsi[errorsi == 0] = 500

        sum_agesu = 0
        sum_agesi = 0

    
        for e in range(0, len(agesu)):
            ynh = (1/(errorsu[e]*np.sqrt(2*np.pi)))*(np.exp(-0.5*(((xi-agesu[e])/errorsu[e])**2)))
            if norm_camels:
                ynh /= np.max(ynh)

            sum_agesu = sum_agesu + ynh

        for o in range(0, len(agesi)):
            ysh = (1/(errorsi[o]*np.sqrt(2*np.pi)))*(np.exp(-0.5*(((xi-agesi[o])/errorsi[o])**2)))
            if norm_camels:
                ysh /= np.max(ysh)

            sum_agesi = sum_agesi + ysh
        
        sum_agesu_list.append(sum_agesu)
        sum_agesi_list.append(sum_agesi)

    vlines = np.arange(tmin,tmax,500)
    pl.ion()
    pl.figure()
    pl.xlim(tmin,tmax)
    pl.xlabel('age [a]')
    pl.ylabel('Relative Probability')
    pl.yticks([])
    pl.title(f'Pseudo time-series summed pdfs RANDOMLY ASIGNED moraine ages NH (blue), SH (red): nsims = {len(rand)}')
    pl.ylim(0,(1.1 * np.max(sum_agesu)))
    pl.vlines(vlines, ymin=0, ymax=(1.1 * np.max((sum_agesu, sum_agesi))),linestyles='dashed')
    for u in range(0, len(sum_agesu_list)):
        pl.plot(xi,sum_agesu_list[u], 'b')
    for i in range(0, len(sum_agesi_list)):
        pl.plot(xi, sum_agesi_list[i], 'r')

    return sum_agesu_list, sum_agesi_list



def obsCorrCoef(y1,y2):
    ynh = y1
    ysh = y2
    n = np.arange(10,1000,10)
    offset_listu = []
    offset_listi = []
    CorrCoef_listu = []
    CorrCoef_listi = []

    for u in range(0, len(n)):
        nu = n[u]
        offsetu = -2*nu
        offset_listu.append(offsetu)

        ynhu = np.concatenate((ynh[nu:], np.zeros(nu)))
        yshu = np.concatenate((np.zeros(nu), ysh[:-nu]))

        ynhuf = ynhu[nu:-nu]
        yshuf = yshu[nu:-nu]

        CorrMatu = np.corrcoef(ynhuf, yshuf)
        CorrCoefu = CorrMatu[0,1]
        CorrCoef_listu.append(CorrCoefu)
    
    for i in range(0, len(n)):
        ni = n[i]
        offseti = 2*ni
        offset_listi.append(offseti)

        ynhi = np.concatenate((np.zeros(ni), ynh[:-ni]))
        yshi = np.concatenate((ysh[ni:], np.zeros(ni)))

        ynhif = ynhi[ni:-ni]
        yshif = yshi[ni:-ni]

        CorrMati = np.corrcoef(ynhif, yshif)
        CorrCoefi = CorrMati[0,1]
        CorrCoef_listi.append(CorrCoefi)

    offset_list = offset_listu[::-1] + offset_listi
    CorrCoef_list = CorrCoef_listu[::-1] + CorrCoef_listi

    opt_offset = offset_list[np.argmax(CorrCoef_list)]

    pl.figure()
    pl.plot(offset_list,CorrCoef_list)
    pl.vlines(offset_list[np.argmax(CorrCoef_list)], -1, 1, colors='g')
    pl.ylim(-1, 1)
    pl.xlabel('Offset [years]')
    pl.ylabel('Correlation Coefficient')
    pl.title(f'Time-series lateral offset CorrCoef optimization: Optimized offset = {offset_list[np.argmax(CorrCoef_list)]}, Max CorrCoef = {round(np.max(CorrCoef_list), 4)}')
    
    
    print(f'Optimized offset = {offset_list[np.argmax(CorrCoef_list)]} with a CorrCoef = {round(np.max(CorrCoef_list), 4)}')
    print('-----------------------')

    return ynh, ysh, opt_offset

def obsCorrCoefRAND(list1,list2, obs_offset1, obs_offset2):
    ynh = list1
    ysh = list2
    n = np.arange(10,1000,10)
    
    master_offlist = []
    master_cclist = []
    
    for h in range(0, len(ynh)):    
        y1 = ynh[h]
        y2 = ysh[h]
        
        offset_listu = []
        offset_listi = []
        CorrCoef_listu = []
        CorrCoef_listi = []

        

        for u in range(0, len(n)):
            nu = n[u]
            offsetu = -2*nu
            offset_listu.append(offsetu)

            ynhu = np.concatenate((y1[nu:], np.zeros(nu)))
            yshu = np.concatenate((np.zeros(nu), y2[:-nu]))

            ynhuf = ynhu[nu:-nu]
            yshuf = yshu[nu:-nu]

            CorrMatu = np.corrcoef(ynhuf, yshuf)
            CorrCoefu = CorrMatu[0,1]
            CorrCoef_listu.append(CorrCoefu)
        
        for i in range(0, len(n)):
            ni = n[i]
            offseti = 2*ni
            offset_listi.append(offseti)

            ynhi = np.concatenate((np.zeros(ni), y1[:-ni]))
            yshi = np.concatenate((y2[ni:], np.zeros(ni)))

            ynhif = ynhi[ni:-ni]
            yshif = yshi[ni:-ni]

            CorrMati = np.corrcoef(ynhif, yshif)
            CorrCoefi = CorrMati[0,1]
            CorrCoef_listi.append(CorrCoefi)

        offset_list = offset_listu[::-1] + offset_listi
        CorrCoef_list = CorrCoef_listu[::-1] + CorrCoef_listi

        
        max_cc = np.max(CorrCoef_list)
        opt_offset = offset_list[np.argmax(CorrCoef_list)]

        master_offlist.append(opt_offset)
        master_cclist.append(max_cc)

    mn_offset = np.mean(master_offlist)
    err_offset = np.std(master_offlist)
    x_fill = np.arange((mn_offset-err_offset), (mn_offset + err_offset), 1)
    moff_array = np.array(master_offlist)

    xi = np.arange(np.min(moff_array) * 1.1,np.max(moff_array) * 1.1,1)
    yi = (1/(err_offset*np.sqrt(2*np.pi)))*(np.exp(-0.5*(((xi-mn_offset)/err_offset)**2)))

    
    counts, bins, _ = pl.hist(moff_array, bins=20, edgecolor='black', linewidth=1)
    
    scalar = np.max(counts) / np.max(yi)

    #hist, bins = np.histogram(moff_array, bins=20, density=True)
    #hist_max = np.max(hist)
    #yi_max = np.max(yi)
    #scalar = yi_max / hist_max

    #portion of code that calculates the z-score and p-value for the observed offset.
    #if p-value is below 0.05, that means that less than 5% of the synthetic data
    #is outside the observed offset and it is very unlikely that the the observed
    #offset would be commonly replicated in the synthetic data.
    obs_offset_filter = obs_offset1
    obs_offset_ALL = obs_offset2

    z_score1 = (obs_offset_filter - mn_offset) / err_offset
    p_val1 = 2 * (1 - stats.norm.cdf(abs(z_score1)))

    z_score2 = (obs_offset_ALL - mn_offset) / err_offset
    p_val2 = 2 * (1 - stats.norm.cdf(abs(z_score2)))

    print(f'p-value for filtered offset = {p_val1}')
    print(f'p-value for all offset = {p_val2}')


    pl.ion()
    pl.figure()
    #pl.plot(xi,(yi * scalar),'k')
    pl.hist(moff_array, bins=20, edgecolor='black', linewidth=1)
    pl.fill_between(x_fill, 0, (np.max(yi * scalar)*1.1), alpha = 0.5)
    #pl.ylim(0,(np.max(yi)*1.1))
    #pl.yticks([])
    pl.xlim(np.min(moff_array) * 1.1,np.max(moff_array) * 1.1)
    pl.xlabel('offset [yrs]')
    pl.ylabel('Count')
    pl.title(f'Histogram of monte carlo-generated offsets, mean = {round(mn_offset, 4)}, std = {round(err_offset, 4)}')

    print(f"Randomly assigned ages offset = {np.mean(master_offlist)} +/ {round(np.std(master_offlist), 4)} with a correlation coefficient = {round(np.mean(master_cclist), 4)} +/ {round(np.std(master_cclist), 4)}")


    return ynh, ysh


###

# All the functions are built above, and then commands below do the analyses
# and make the plots in an organized way

###


#STEP 1: Get results from ICE-D Queries
#first one is for the 'filtered dataset'
#second one is for the complete, non filtered dataset
query_result = bps_query()
query_resultALL = bps_queryALL()

#STEP 2: Go thru analysis of filtering 'filtered dataset'
group_result = grouping(query_result)
filtered_result = Chi2_filtering(group_result)

#STEP 3a: Organize data from ICE-D and plot 'filtered dataset' as summed pdf
y1, y2 = data_organize(filtered_result)
#STEP 3b: Determine the optimal offset to get the best correlation coefficient for the 'filtered dataset'
ynh, ysh, ofst1 = obsCorrCoef(y1,y2)

#STEP 4a: Organize data from ICE-D and plot 'full dataset' as summed pdf
y3, y4 = data_organizeALL(query_resultALL)
#STEP 4b: Determine the optimal offset to get the best correlation coefficient for the 'filtered dataset'
ynhALL, yshALL, ofst2 = obsCorrCoef(y3,y4)

#STEP 5a: Organize data from ICE-D
#run monte carlo simulation to randomly assign ages to two new datasets iteratively
#plot ensemble of monte carlo simulations as summed pdfs
list1, list2 = data_organizeRAND(query_resultALL)

#STEP 5b: For each simulation, determine the optimal offset to get the best correlation coefficient
#plot a normal pdf using the average and standard deviation of the simulated offset to achieve max Corr Coef
yrand1, yrand2 = obsCorrCoefRAND(list1, list2, ofst1, ofst2)