import urllib.request
import os
from urllib.error import URLError,HTTPError
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
	
	def check_path(self,dir):
		'''
		This fucntion check if this directory exists, otherwise it will create one
		- dir the directory to check or created.
		- return True if the path has been created successfully
		'''
		if not os.path.exists(dir):
			try:
				os.makedirs(dir) # check if the directory exists, create one if not
				return True
			except OSError as e: # Trap the OS error and show the embedded error code
				if e.errno != errno.EEXIST:
					raise
    				
	def download_files(self,url,dir,file = None):
		'''
		Download file from given url to given directory with all handled errors
		- url The string that represents correct url that has the file
		- dir a local directory that exists 
		'''
		#if file is not None:
			
		try:
			urllib.request.urlretrieve(url,dir)
		except URLError as e:
			print("Invalid URL: " + e.reason)
		except HTTPError as e:
			print("HTTP ERROR:" + e.code)