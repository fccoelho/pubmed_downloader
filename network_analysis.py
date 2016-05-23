
import networkx as nx
import pymongo
from gensim import corpora, models, similarities
from nltk.tokenize import WordPunctTokenizer
from nltk.corpus import stopwords
from string import punctuation

connection = pymongo.MongoClient()
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

def article_generator(artcol, tokenize=True):
    cursor_art = connection.pubmed[artcol].find()
    for art in cursor_art:
        try:
            abs = art["MedlineCitation"]["Article"]['Abstract']['AbstractText']
        except KeyError:
            abs = ''
        if tokenize:
            yield [token for token in tokenizer.tokenize(str(abs)) if token not in sw]
        else:
            yield abs

def create_article_dictionary(collection='zika'):
    artcol = "articles" if collection == 'zika' else collection
    articles = article_generator(artcol)
    dictionary = corpora.Dictionary(articles)
    dictionary.save('Dicionario_{}.dict'.format(collection))
    return dictionary

def create_corpus(collection='zika'):
    artcol = "articles" if collection == 'zika' else collection
    dictionary = create_article_dictionary(collection)
    # print(dictionary)
    articles = article_generator(artcol)
    corpus = [dictionary.doc2bow(text) for text in articles]
    corpora.MmCorpus.serialize('corpus_{}'.format(collection), corpus)
    return corpus, dictionary

def LSI_topics(corpus, dictionary):
    tfidf = models.TfidfModel(corpus)
    corpus_tfidf = tfidf[corpus]
    #print(corpus_tfidf)
    lsi = models.LsiModel(corpus_tfidf, id2word=dictionary, num_topics=10)
    return lsi

def LDA_topics(corpus, dictionary, num_topics):
    lda_model = models.ldamodel.LdaModel(corpus, id2word=dicionary, num_topics=num_topics, passes=10)
    return lda_model




if __name__ == "__main__":
    # dump_to_csv()
    # dump_to_csv('mers')
    col = 'zika'
    c, d = create_corpus(col)
    lsi = LSI_topics(c, d)
    lsi.save('lda_model_{}'.format(col))
    lda = LDA_topics(c, d, 30)
    lda.save('lda_model_{}'.format(col))
    print(lsi.show_topics(10))


