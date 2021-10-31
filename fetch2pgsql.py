#!/bin/env python3
"""
Fetches records from pubmed database stores them into a postgresql table and enables full-text indexing on it
"""


from Bio import Entrez
from urllib.error import HTTPError, URLError
from sqlalchemy import create_engine
import time
from tqdm import tqdm
from getpass import getpass
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
from datetime import date


pw = getpass("Enter database password: ")

engine = create_engine(f"postgresql://postgres:{pw}@localhost/pubmed")


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
months = {'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6,'Jul':7, 'Aug':8, 'Sep':9, 'Oct':10, 'Nov':11, 'Dec':12}
class SearchAndCapture:
    def __init__(self, email, search_term, collection='articles'):
        self.search_term =search_term
        Entrez.email = email
        self.articles_table = "articles"
        self.citation_table = "citations"
        with engine.connect() as connection:
            connection.execute("create table IF NOT EXISTS articles(pmid integer unique, title varchar(1024), abstract text, journal varchar(512), pubdate date);")
            connection.execute("create table IF NOT EXISTS citations(pmid integer unique, cited_by text);")
            connection.execute(f"create unique index IF NOT EXISTS pmid_idx on {self.articles_table} (pmid)")
            connection.execute(f"create unique index IF NOT EXISTS pmid_idx on {self.citation_table} (pmid)")
        # self.collection = connection.pubmed[collection]
        # self.collection.create_index([('$**', pymongo.TEXT)], default_language='english')
        # if collection == 'articles':
        #     self.citation_colection = connection.pubmed.citations
        # else:
        #     self.citation_colection = connection.pubmed["citations_{}".format(collection)]

    def _fetch(self, pmid):
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode='xml')
        records = Entrez.read(handle)
        return records.get('PubmedArticle',[])

    def _get_old_ids(self):
        with engine.connect() as connection:
            res = connection.execute(f"select pmid from {self.articles_table};")
            oldids = res.fetchall()
        
        return [i[0] for i in oldids]

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

    def create_FT_index(self):
        with engine.connect() as connection:
            query = f"ALTER TABLE articles ADD COLUMN textsearchable_index_col tsvector\
               GENERATED ALWAYS AS (to_tsvector('english', coalesce(title, '') || ' ' || coalesce(abstract, ''))) STORED;"
            query2 = "CREATE INDEX textsearch_idx ON articles USING GIN (textsearchable_index_col);"
            connection.execute(query)

    def update_citations_concurrently(self):
        print("Updating citations...")
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
                        with engine.connect() as connection:
                            cits = '|'.join(cits)
                            connection.execute(f"insert into {self.citation_table} VALUES(%s, %s) \
                            on conflict (pmid) do update set cited_by=%s;",(pmid, cits,cits))
                except Exception as e:
                    print(f'{pmid} generated an exception: {e}')

    def update_citations(self):
        print("Updating citations...")
        ids = self._get_old_ids()
        for i in ids:
            cits = self._get_citations(i)
            if cits is None:
                continue
            with engine.connect() as connection:
                            cits = '|'.join(cits)
                            connection.execute(f"insert into {self.citation_table} VALUES(%s, %s) \
                            on conflict (pmid) do update set cited_by=%s;",(pmid, cits,cits))

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
        # print (old_ids)
        # return
        print("Fetching results for {}".format(self.search_term))
        handle = Entrez.esearch(db="pubmed", retmax=100000, term=self.search_term)
        response = Entrez.read(handle=handle)
        print("Found {} items".format(len(response['IdList'])))
        new_ids = {}
        for pmid in tqdm(response['IdList']):
            if int(pmid) in old_ids:
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
            with engine.connect() as connection:
                article = art["MedlineCitation"]["Article"]
                try:
                    abstract = article['Abstract']['AbstractText'][0].replace("'","''")
                except (AttributeError, KeyError) as e:
                    print(e)
                    abstract = ""
                title = article["ArticleTitle"].replace("'","''")
                journal = article["Journal"]["Title"]
                year = int(article["Journal"]["JournalIssue"]["PubDate"].get('Year', '2016'))
                mo = article["Journal"]["JournalIssue"]["PubDate"].get('Month', 'Jan')
                try: 
                    month = months[mo]
                except KeyError as e:
                    # print(e)
                    month = int(mo[1] if mo.startswith('0') else mo)
                day = int(article["Journal"]["JournalIssue"]["PubDate"].get('Day', '1'))
                data = date(year,month,day)
                query = f"insert into {self.articles_table} VALUES(%s, %s, %s, %s, %s) \
                on conflict (pmid) do update set title=%s, abstract=%s, journal=%s, pubdate=%s;"
                connection.execute(query, (int(pmid),title, abstract, journal, data.isoformat(), title, abstract, journal, data.isoformat()))
            new_ids[pmid] = pmid
        return new_ids





if __name__ == "__main__":
    S = SearchAndCapture('fccoelho@gmail.com', '((zika microcephaly) NOT zika[author])')
    S.update_multiple_searches(zika_query_strings)
    S.update_citations_concurrently()
    S.create_FT_index()
    # for s in [MERS_query_strings, Mayaro_query_strings, Oropouche_query_strings]:
    #     T = SearchAndCapture('fccoelho@gmail.com', s[0], s[0].lower())
    #     T.update_multiple_searches(s)
    #     T.update_citations_concurrently()
