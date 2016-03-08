__author__ = 'fccoelho'


import scholarly
import pymongo
import time
import json

conn = pymongo.MongoClient()


search_query = scholarly.search_pubs_query('zika zikv -author:zika')

def continuous_fetch():
    while True:
        doc = {}
        try:
            art = next(search_query)
            if not art._filled:
                art.fill()
            doc['bib'] = art.bib
            doc['citedby'] = art.citedby
            doc['id_scholarcitedby'] = art.id_scholarcitedby
            doc['url_scholarbib'] = art.url_scholarbib
            doc['source'] = art.source
            conn.scholar.articles.update_one({'url_scholarbib': art.url_scholarbib}, {"$setOnInsert": doc}, upsert=True)
        except Exception as e:
            print(e, art)
            break
        time.sleep(0.5)

if __name__ == "__main__":
    continuous_fetch()
