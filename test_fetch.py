#-*- coding:utf-8 -*-
__author__ = 'fccoelho'

import unittest
from fetch import SearchAndCapture
from bson.objectid import ObjectId


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.S = SearchAndCapture('fccoelho@gmail.com', 'zika virus')

    def test_fetch_one_article(self):
        resp = self.S._fetch('26923117')
        self.assertIsInstance(resp, list)
        self.assertEqual(len(resp), 1)

    def test_get_old_ids(self):
        res = self.S._get_old_ids()
        self.assertIsInstance(res[0], str)


    def test_update_database(self):
        ids = self.S.update()
        self.assertIsInstance(ids, dict)
        self.assertGreater(len(ids), 0)
        obj = list(ids.values())[0]
        self.assertIsInstance(obj, ObjectId)



if __name__ == '__main__':
    unittest.main()
