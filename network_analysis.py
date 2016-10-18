
import networkx as nx
import pymongo
from gensim import corpora, models, similarities
from nltk.tokenize import WordPunctTokenizer
from nltk.corpus import stopwords
from string import punctuation

connection = pymongo.MongoClient('172.16.4.51')
tokenizer = WordPunctTokenizer()
sw = stopwords.words('english') + list(punctuation)

def dump_to_csv(collection=""):
    artcol = "articles" if not collection else collection
    citcol = "citations" if not collection else "citations_"+collection
    cursor_cit = connection.pubmed[citcol].find()
    cursor_art = connection.pubmed[artcol].find()
    f = open("{}edge_list.csv".format(collection), 'w')
    g = open("{}node_list.csv".format(collection),'w')
    g.write("PMID|title|journal|year|month|day\n")
    try:
        for cit in cursor_cit:
            # if cit['citedby']:
            for id_citer in cit['citedby']:
                f.write("{},{}\n".format(id_citer, cit['PMID']))
        for art in cursor_art:
            article = art["MedlineCitation"]["Article"]
            g.write("{}|{}|{}|{}|{}|{}\n".format(art["MedlineCitation"]["PMID"],
                                                 article["ArticleTitle"],
                                                 article["Journal"]["Title"],
                                                 article["Journal"]["JournalIssue"]["PubDate"].get('Year', '2016'),
                                                 article["Journal"]["JournalIssue"]["PubDate"].get('Month', 'NA'),
                                                 article["Journal"]["JournalIssue"]["PubDate"].get('Day', 'NA'),
                                                 ))
    finally:
        f.close()
        g.close()

def article_generator(artcol, filter={}, tokenize=True):
    cursor_art = connection.pubmed[artcol].find(filter)
    for art in cursor_art:
        try:
            abs = art["MedlineCitation"]["Article"]['Abstract']['AbstractText']
        except KeyError:
            abs = ''
        if tokenize:
            yield [token for token in tokenizer.tokenize(str(abs)) if token not in sw]
        else:
            yield abs

def create_article_dictionary(collection='zika', filter=None):
    if filter is None:
        filter = {}
    artcol = "articles" if collection == 'zika' else collection
    articles = article_generator(artcol, filter=filter)
    dictionary = corpora.Dictionary(articles)
    dictionary.save('Dicionario_{}.dict'.format(collection))
    return dictionary

def create_corpus(collection='zika', filter=None):
    if filter is None:
        filter = {}
    artcol = "articles" if collection == 'zika' else collection
    dictionary = create_article_dictionary(collection, filter=filter)
    stopw_ids = map(dictionary.token2id.get, sw)
    dictionary.filter_tokens(stopw_ids)
    dictionary.compactify()
    dictionary.filter_extremes(no_below=5, no_above=0.5, keep_n=None)
    dictionary.compactify()
    # print(dictionary)

    articles = article_generator(artcol, filter=filter)
    corpus = [dictionary.doc2bow(text) for text in articles]
    corpora.MmCorpus.serialize('corpus_{}'.format(collection), corpus)
    return corpus, dictionary

def get_top_topics_by_year(lda, collection, year):
    year = str(year)
    corpus, dic = create_corpus(collection=collection, filter={"MedlineCitation.DateCreated.Year": year})
    corpus_lda = lda[corpus]
    return [sorted(doc, key=lambda item: -item[1]) for doc in corpus_lda]

def LSI_topics(corpus, dictionary):
    tfidf = models.TfidfModel(corpus)
    corpus_tfidf = tfidf[corpus]
    #print(corpus_tfidf)
    lsi = models.LsiModel(corpus_tfidf, id2word=dictionary, num_topics=10)
    return lsi

def LDA_topics(corpus, dictionary, num_topics):
    lda_model = models.ldamodel.LdaModel(corpus, id2word=dictionary, num_topics=num_topics, passes=10)
    return lda_model


if __name__ == "__main__":
    dump_to_csv()
    dump_to_csv('mers')
    for col in ['zika', 'mers']:
        print("Calculating {} lsi model".format(col))
        c, d = create_corpus(col)
        lsi = LSI_topics(c, d)
        print("Saving {} lsi model".format(col))
        lsi.save('lda_model_{}'.format(col))
        print("Calculating {} lda model".format(col))
        lda = LDA_topics(c, d, 30)
        print("Saving {} lda model".format(col))
        lda.save('lda_model_{}'.format(col))
    # print(lsi.show_topics(10))


