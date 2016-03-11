__author__ = 'fccoelho'


import scholarly
import pymongo
import time
import json

conn = pymongo.MongoClient()


search_query = scholarly.search_pubs_query('zika zikv -author:zika')

def continuous_fetch():
    downloaded = [a['url_scholarbib'] for a in conn.scholar.articles.find({}, {'url_scholarbib': 1})]
    while True:
        doc = {}
        # try:
        art = next(search_query)
        if art.url_scholarbib in downloaded:
            continue
        if not art._filled:
            art.fill()
        doc['bib'] = art.bib
        try:
            doc['citedby'] = art.citedby
        except AttributeError:
            doc['citedby'] = 0
        try:
            doc['id_scholarcitedby'] = art.id_scholarcitedby
        except AttributeError:
            doc['id_scholarcitedby'] = ""
        doc['url_scholarbib'] = art.url_scholarbib
        doc['source'] = art.source
        conn.scholar.articles.update_one({'url_scholarbib': art.url_scholarbib}, {"$setOnInsert": doc}, upsert=True)
        print("inserted {}".format(art.url_scholarbib))
        # except Exception as e:
        #     print("==>", e)
        #     print(art)
        #     break
        time.sleep(0.1)

if __name__ == "__main__":
    continuous_fetch()
