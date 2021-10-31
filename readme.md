# Pubmed Downloader

This code helps you to start a thematic collection of articles by automatically downloading them from PubMed.

It uses Pubmed official API, so there is no scraping involved and its completely in accordance with its usage policy.

## Downloading

To start using it, follow the steps below.

```python
from fetch import SearchAndCapture
S = SearchAndCapture('your@email.com', '((zika microcephaly) NOT zika[author])')
S.update_multiple_searches()
S.update_citations_concurrently()
```

That's all! You may want to run it once a day to include recently published articles. This will also update a separate collection with the citation network for all the `PMID`s in the article collection

You need to have a local mongodb server running.

## Downloading to postgresql
If you don't want to work with Mongodb,

You can create your corpus of articles instead on PostgreSQL. 

For that you can customize your query strings on the `fetch2psql.py` script and run it on the command line.

```bash
$ ./fetch2pgsql.py
```

This script also enables full text indexing on the corpus.

## Export network data
You can export csv files with data structured for network analysis:

```
$ python3 network_analysis.py
```
This script only works for the mongodb databases. I should not be difficult to adapt it to work with corpora stored in Postgres.