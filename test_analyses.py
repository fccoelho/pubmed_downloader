import unittest
from gensim import models
from network_analysis import create_corpus, get_top_topics_by_year, LDA_topics


class TestAnalyses(unittest.TestCase):
    def test_topics_by_year(self):
        lda = models.LdaModel.load('lda_model_zika')
        c,d = corpus, dic = create_corpus('zika')
        dt = get_top_topics_by_year(lda, d, 'zika', 2016)

        self.assertIsInstance(dt, list)
    def test_lda_topics(self):
        c, d = corpus, dic = create_corpus('zika')
        lda = LDA_topics(c, d, 15)




if __name__ == '__main__':
    unittest.main()
