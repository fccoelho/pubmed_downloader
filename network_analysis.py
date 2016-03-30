
import networkx as nx
import pymongo

connection = pymongo.MongoClient()

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

if __name__ == "__main__":
    dump_to_csv()
    dump_to_csv('mers')
