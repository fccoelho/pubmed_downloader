
u"""
Created on 03/03/16
by fccoelho
license: GPL V3 or Later
"""

__docformat__ = 'restructuredtext en'
from Bio import Entrez
import pymongo

connection = pymongo.MongoClient()


class SearchAndCapture:
    def __init__(self, email, search_term):
        self.search_term =search_term
        Entrez.email = email
        self.collection = connection.pubmed.articles

    def update(self):
        handle = Entrez.esearch("pubmed", term=self.search_term)
        response = Entrez.read(handle=handle)


