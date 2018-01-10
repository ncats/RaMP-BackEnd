import urllib.request
import os
import time
import datetime 
class MetabolomicsData():
	'''
	This class is a super class for hmdb, kegg, wiki, Reactome database,
	which has the general functions defined for all other classes
	'''
	def __init__(self):
		self.day = datetime.datetime.today()
		#self.day = "Hello World"
	def __str__(self):
		return self.date
		
	def download_files(self,url,dir):
		urllib.request.urlretrieve(url,dir)
