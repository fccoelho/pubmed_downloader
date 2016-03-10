
import networkx as nx
import pymongo

connection = pymongo.MongoClient()

def dump_to_csv():
    cursor_cit = connection.pubmed.citations.find()
    cursor_art = connection.pubmed.articles.find()
    f = open("edge_list.csv", 'w')
    g = open("node_list",'w')
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
