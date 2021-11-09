
#!/bin/env python3

import networkx as nx
from nltk import corpus
from sqlalchemy import create_engine
from gensim import corpora, models, similarities
from gensim.models import Word2Vec, KeyedVectors, FastText, Doc2Vec
from gensim.models.ldamulticore import LdaMulticore
from gensim.utils import save_as_line_sentence
from gensim.models.word2vec import LineSentence
from gensim.models.doc2vec import TaggedLineDocument
from nltk.tokenize import WordPunctTokenizer, PunktSentenceTokenizer
from nltk.corpus import stopwords
from string import punctuation
from getpass import getpass
import pandas as pd


pw = getpass("Enter database password: ")
engine = create_engine(f"postgresql://postgres:{pw}@localhost/pubmed",pool_size=20, max_overflow=0)
tokenizer = WordPunctTokenizer().tokenize
sent_tokenizer = PunktSentenceTokenizer().tokenize
sw = stopwords.words('english') + list(punctuation)

def dump_to_csv(table=""):
    arttbl = "articles" if not table else table
    cittbl = "citations" if table=="articles" else "citations_"+table
    with engine.connect() as conn:
        dfa = pd.read_sql_table(arttbl, conn, chunksize=1000)
        dfc = pd.read_sql_table(cittbl, conn, chunksize=1000)
    header=True
    for chunk in dfa:
        chunk.to_csv(f"{arttbl}_node_list.csv", header=header, mode='a')
        header=False
    header=True
    for chunk in dfc:
        chunk.to_csv(f"{cittbl}_edge_list.csv", header=header, mode='a')
        header=False



def abstract_generator(artcol, sql_query=None):
    if sql_query is None:
        sql_query = f"select abstract from {artcol};"
    with engine.connect() as conn:
        res = conn.execute(sql_query)
        for abs in res.fetchall():
            if abs[0].strip() == "":
                continue
            else:
                yield abs[0]

def sentence_generator(texts):
    '''Yields documents as a list of tokenized sentences'''
    for t in texts:
        yield sent_tokenizer(t)


def token_generator(texts):
    '''Yields documents as a list of tokens'''
    for t in texts:
        yield tokenizer(t)


def dump_to_line_sequence_file(table):
    docs = abstract_generator(table)
    tokens = token_generator(docs)
    save_as_line_sentence(tokens, 'sentences.txt')

def calculate_w2v(corpus_file):
    model_sent = Word2Vec(sentences=LineSentence(corpus_file), workers=8)
    model_sent.save('abstracs.w2v')

def calculate_d2v(corpus_file):
    model_sent3 = Doc2Vec(documents=TaggedLineDocument(corpus_file), workers=8)
    model_sent3.save('abstracts.d2v')

def create_article_dictionary(table='articles'):
    articles = abstract_generator(table)
    tokens = token_generator(articles)
    dictionary = corpora.Dictionary(tokens)
    dictionary.save(f'Dicionario_{table}.dict')
    return dictionary

def create_corpus(table='articles'):
    dictionary = create_article_dictionary(table)
    stopw_ids = map(dictionary.token2id.get, sw)
    dictionary.filter_tokens(stopw_ids)
    dictionary.compactify()
    dictionary.filter_extremes(no_below=5, no_above=0.5, keep_n=None)
    dictionary.compactify()
    # print(dictionary)

    articles = abstract_generator(table)
    corpus = [dictionary.doc2bow(tokenizer(text)) for text in articles]
    corpora.MmCorpus.serialize('corpus_{}'.format(table), corpus)
    return corpus, dictionary

def get_top_topics_by_year(lda, full_dict, table, year):
    year = str(year)
    articles = abstract_generator(table, sql_query=f"select abstract from {table} where date_part('year',pubdate)={year}")
    corpus = [full_dict.doc2bow(tokenizer(text)) for text in articles]
    year_lda = lda[corpus]
    print(lda.top_topics(corpus))
    return lda.top_topics(corpus) #sorted([doc for doc in corpus_lda], key=lambda item: item[1], reverse=True)

def LSI_topics(corpus, dictionary):
    tfidf = models.TfidfModel(corpus)
    corpus_tfidf = tfidf[corpus]
    #print(corpus_tfidf)
    lsi = models.LsiModel(corpus_tfidf, id2word=dictionary, num_topics=10)
    return lsi

def LDA_topics(corpus, dictionary, num_topics):
    lda_model = LdaMulticore(corpus, id2word=dictionary, num_topics=num_topics, passes=10)
    return lda_model


if __name__ == "__main__":
    # dump_to_csv()
    dump_to_csv('articles')
    dump_to_line_sequence_file('articles')
    calculate_w2v('sentences.txt')
    # for col in ['articles']:
    #     print("Calculating {} lsi model".format(col))
    #     c, d = create_corpus(col)
    #     lsi = LSI_topics(c, d)
    #     print("Saving {} lsi model".format(col))
    #     lsi.save('lda_model_{}'.format(col))
    #     print("Calculating {} lda model".format(col))
    #     lda = LDA_topics(c, d, 30)
    #     print("Saving {} lda model".format(col))
    #     lda.save('lda_model_{}'.format(col))
    #     doc_topics = get_top_topics_by_year(lda, d, 'articles', 2016)
    # print(lsi.show_topics(10))


