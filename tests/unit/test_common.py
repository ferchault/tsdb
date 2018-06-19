#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import sys
import os

# load tsdb
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'src'))
import tsdb

class TestCommon(unittest.TestCase):
	def test_get_elements_from_label_fails_on_unkown(self):
		with self.assertRaises(ValueError):
			tsdb.get_elements_from_label('foobar_A_A_A_A_R_R')

	def test_get_elements_from_label_fails_on_wrong_length(self):
		with self.assertRaises(ValueError):
			tsdb.get_elements_from_label('foobar_A_A_A_A_R_R_X_X_X')

	def test_get_elements_from_label_sn2(self):
		ret = tsdb.get_elements_from_label('sn2_C_D_E_D_C_B')
		expected = set('C H Br F C C N C H H H N H H C H H H'.split())
		self.assertEqual(ret, expected)

	def test_get_elements_from_label_e2(self):
		ret = tsdb.get_elements_from_label('e2ts_D_C_A_D_C_A')
		expected = set('C C Br H H H C H H H C N C H H H '.split())
		self.assertEqual(ret, expected)						

