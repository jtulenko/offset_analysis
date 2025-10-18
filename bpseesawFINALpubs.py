import pymysql
import numpy as np
import pandas
import warnings
import bibtexparser
from tabulate import tabulate

import db_info

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

def bps_queryALL():
    dbc = reader_connect_to_db()
    dbcursor = dbc.cursor()

    query = f"""SELECT DISTINCT base_publication.bibtex_record, base_publication.doi
        FROM base_sample
        JOIN base_site ON base_sample.site_id = base_site.id
        JOIN base_calculatedage ON base_calculatedage.sample_id = base_sample.id
        JOIN base_application_sites ON base_application_sites.site_id = base_site.id
        LEFT JOIN base_samplepublicationsmatch ON base_sample.id = base_samplepublicationsmatch.sample_id
        LEFT JOIN base_publication ON base_samplepublicationsmatch.publication_id = base_publication.id
        WHERE base_site.what LIKE "%moraine%"
        AND ({region} OR base_sample.lat_DD < -30)
        AND base_application_sites.application_id = 2
        AND (base_calculatedage.t_{age_scale} > {tmin - 1000} AND base_calculatedage.t_{age_scale} < {tmax + 1000})"""
    dbcursor.execute(query)
    result = dbcursor.fetchall()
    
    dbc.close()

    result_list = np.array(result)

    return result_list

def publication_list(result_list):
    result_list = result_list

    bib_list = result_list[:,0]

    authors=[]
    years=[]
    titles=[]
    journals=[]
    volumes=[]
    numbers=[]
    pages=[]
    dois = result_list[:,1]

    for record in bib_list:
        try:
            bib = bibtexparser.loads(record)
            entry = bib.entries[0]
            authors.append(entry.get("author", ""))
            years.append(entry.get("year", ""))
            titles.append(entry.get("title", ""))
            journals.append(entry.get("journal", ""))
            volumes.append(entry.get("volume", ""))
            numbers.append(entry.get("number", ""))
            pages.append(entry.get("pages", ""))

        except Exception as e:
            print(f"Error: {e}")
            years.append(entry.get(""))
            titles.append(entry.get(""))
            journals.append(entry.get(""))
            volumes.append(entry.get(""))
            numbers.append(entry.get(""))
            pages.append(entry.get(""))
    
    formatted_list = [f"{a}. ({b}). {c}. {d}, {e}({f}), {g}. {h})" for a,b,c,d,e,f,g,h in zip(authors, years, titles, journals, volumes, numbers, pages, dois)]
    formatted_list.sort()
    #print(formatted_list)

    with open("/home/jtulenko/sw/analyses/bpseesaw_AGU23/citations.txt", "w", encoding="utf-8") as f:
        for entry in formatted_list:
            f.write(entry + "\n")

    # boneyard of code that didn't work
    # df = pandas.DataFrame({
    #     "authors": authors,
    #     "years": years,
    #     "titles": titles,
    #     "journals": journals,
    #     "volumes":volumes,
    #     "numbers": numbers,
    #     "pages": pages,
    #     "DOIs": dois
    # })

    #table = list(zip(authors, years, titles, journals, volumes, numbers, pages, dois))
    #print(tabulate(table, headers=["authors", "years", "titles", "journals", "volumes", "numbers", "pages", "DOIs"]))

    return result_list

query_resultALL = bps_queryALL()
test = publication_list(query_resultALL)