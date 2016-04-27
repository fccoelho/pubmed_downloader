# Pubmed Downloader

This code helps you to start a thematic collection of articles by automatically downloading them from PubMed.

It uses Pubmed official API, so there is no scraping involved and its completely in accordance with its usage policy.

To start using it, follow the steps below.

```python
S = SearchAndCapture('your@email.com', '((zika microcephaly) NOT zika[author])')
S.update_multiple_searches(zika_query_strings)
S.update_citations_concurrently()
```

That's all! You may want to run update() once a day to include recently published articles.

You need to have a local mongodb server running.
