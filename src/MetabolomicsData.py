import urllib.request
import os
from urllib.error import URLError,HTTPError
import time
import datetime 
from django.core.checks import database
from builtins import str
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
		- file name of the file you want to downloaded. 
			if not None: need to be concatenated with directory {dir}
		'''
		#if file is not None:
		if file is None:	
			try:
				urllib.request.urlretrieve(url,dir)
			except URLError as e:
				print("Invalid URL: " + e.reason)
			except HTTPError as e:
				print("HTTP ERROR:" + e.code)
		elif isinstance(file,str):
			if file not in os.listdir(dir):
				print("{}. Downloading {} ".format(str(len(os.listdir(dir))),file))
				try:
					urllib.request.urlretrieve(url,dir+file)
				except URLError as e:
					print("Invalid URL: " + e.reason)
				except HTTPError as e:
					print("HTTP ERROR:" + e.code)
		else:
			raise TypeError("Expect a string for the file name")	
		
	def write_myself_files(self,database,dir = '../misc/output/'):
		path = dir + database +'/'
		self.check_path(path)
		attrs = vars(self)
		for key in attrs:
			if type(attrs[key]) is dict and len(attrs[key]) > 0:
				with open(path+database+key+".txt",'wb') as f:
					for id, value in attrs[key].items():
					
						if type(value) is not dict and type(value) is not list:
							
							f.write(id.encode('utf-8') +b'\t'
									+ value.encode('utf-8') +b'\n')
						elif type(value) is list:
							for item in value:
								if item is not None and item is not 'NA':
									f.write(id.encode('utf-8') +b'\t'
											+ item.encode('utf-8') +b'\n')
						elif type(value) is dict:
							for source,sourceid in value.items():
								if type(sourceid) is list:
									for each in sourceid:
										if each is not None and each is not 'NA':
											f.write(id.encode('utf-8') +b'\t'
											+ source.encode('utf-8')+b'\t'
											+ each.encode('utf-8') +b'\n')
								else:
									if sourceid is not None and sourceid is not 'NA':
										f.write(id.encode('utf-8') +b'\t'
											+ source.encode('utf-8')+b'\t'
											+ str(sourceid).encode('utf-8') +b'\n')
									
					
				
				
	

		
		