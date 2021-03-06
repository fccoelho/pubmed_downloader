
u"""
Created on 03/03/16
by fccoelho
license: GPL V3 or Later
"""

__docformat__ = 'restructuredtext en'
from Bio import Entrez
from urllib.error import HTTPError, URLError
import pymongo
import time
from concurrent.futures import ThreadPoolExecutor, as_completed


connection = pymongo.MongoClient('172.16.4.51')

zika_query_strings = ['((zika) NOT zika[Author])',
                 '((zika virus NOT zika[Author]))',
                 '((zika microcephaly) NOT zika[author])',
                 '((zika fever NOT zika[Author]))',
                 '(zika NOT zika[Author] ) guillain barre',
                 '((zika epidemic NOT zika[Author]))',
                 '((zika model NOT zika[Author]))',
                 '(zika NOT zika[Author] ) chikungunya',
                 '(zika NOT zika[Author] ) outbreak',
                 '(zika NOT zika[Author] ) dengue',
                 '(zika NOT zika[Author]) culex ',
                 '(zika NOT zika[Author]) mosquito',
                 '(zika NOT zika[Author])virus transmission'
                 ]
MERS_query_strings = [
    'MERS',
    'mers cov',
    'mers coronavirus',
    'mers korea',
    'mers vaccine',
    'mers camel'
]

Mayaro_query_strings = ['mayaro',
                        'mayaro virus',
                        ]
Oropouche_query_strings = ['oropouche',
                           'oropouche virus'
                           ]

class SearchAndCapture:
    def __init__(self, email, search_term, collection='articles'):
        self.search_term =search_term
        Entrez.email = email
        self.collection = connection.pubmed[collection]
        self.collection.create_index([('$**', pymongo.TEXT)], default_language='english')
        if collection == 'articles':
            self.citation_colection = connection.pubmed.citations
        else:
            self.citation_colection = connection.pubmed["citations_{}".format(collection)]

    def _fetch(self, pmid):
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode='xml')
        records = Entrez.read(handle)
        return records.get('PubmedArticle',[])

    def _get_old_ids(self):
        oldids = self.collection.find({}, {"MedlineCitation.PMID": 1})
        oldids = [i['MedlineCitation']['PMID'] for i in oldids]
        return oldids

    def _get_citations(self, pmid):
        try:
            handle = Entrez.elink(dbfrom="pubmed", id=pmid, linkname="pubmed_pubmed_citedin")
        except HTTPError as e:
            print(e)
            return None
        record = Entrez.read(handle)
        handle.close()
        if record[0]['LinkSetDb']:
            citers = record[0]['LinkSetDb'][0]['Link']
        else:
            citers = []
        return [i['Id'] for i in citers]

    def update_citations_concurrently(self):
        ids = self._get_old_ids()
        with ThreadPoolExecutor(max_workers=100) as executor:
            future_cits = {executor.submit(self._get_citations, pmid): pmid for pmid in ids}
            for future in as_completed(future_cits):
                pmid = future_cits[future]
                print("getting citations for {}".format(pmid))
                try:
                    cits = future.result()
                    if cits is None:
                        continue
                    else:
                        self.citation_colection.update_one({"PMID": pmid}, {"$set": {"citedby": cits}}, upsert=True)
                except Exception as e:
                    print('%r generated an exception: %s' % (pmid, e))

    def update_citations(self):
        print("Updating citations...")
        ids = self._get_old_ids()
        for i in ids:
            cits = self._get_citations(i)
            if cits is None:
                continue
            self.citation_colection.update_one({"PMID": i}, {"$set": {"citedby": cits}}, upsert=True)

    def update_multiple_searches(self, queries=None):

        if queries is not None:
            query_strings = queries
        else:
            query_strings = [self.search_term]
        for qs in query_strings:
            self.search_term = qs
            self.update()

    def update(self):
        old_ids = self._get_old_ids()
        print("Fetching results for {}".format(self.search_term))
        handle = Entrez.esearch(db="pubmed", retmax=100000, term=self.search_term)
        response = Entrez.read(handle=handle)
        print("Found {} items".format(len(response['IdList'])))
        new_ids = {}
        for pmid in response['IdList']:
            if pmid in old_ids:
                continue
            try:
                art = self._fetch(pmid)[0]
            except URLError:
                print("Downloading of {} failed. Skipping".format(pmid))
                continue
            except IndexError as e:
                print("Empty record for {}".format(pmid))
                continue
            if "MedlineCitation" not in art:
                continue
            # TODO: Fix the date before insertion
            result = self.collection.update_one({"MedlineCitation.PMID": str(art["MedlineCitation"]["PMID"])}, {"$setOnInsert": art}, upsert=True)

            new_ids[pmid] = result.upserted_id
        return new_ids





if __name__ == "__main__":
    S = SearchAndCapture('fccoelho@gmail.com', '((zika microcephaly) NOT zika[author])')
    S.update_multiple_searches(zika_query_strings)
    S.update_citations_concurrently()
    for s in [MERS_query_strings, Mayaro_query_strings, Oropouche_query_strings]:
        T = SearchAndCapture('fccoelho@gmail.com', s[0], s[0].lower())
        T.update_multiple_searches(s)
        T.update_citations_concurrently()
