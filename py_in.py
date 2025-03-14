import subprocess
import os
import pandas as pd
import copy
import re
import global_constant as gc

def find_path():
	'''
	read the path of QC software that are written in .appspathrc
	'''
	script_path = os.path.abspath(__file__)
	script_dir = os.path.dirname(script_path)
	pathrc = script_dir + "/.appspathrc"
	pathdf = pd.read_csv(pathrc, sep='\t', header=None)
	pathdict = pd.Series(pathdf[1].values, index=pathdf[0]).to_dict()
	return pathdict


def call_multiwfn(file_path, commands):
	multiwfn_path = find_path()['multiwfn']

	fchk = open(file_path, 'rb')

	process = subprocess.Popen([multiwfn_path, file_path], stdin=subprocess.PIPE,
	stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

	stdout, stderr = process.communicate(input=commands)
	fchk.close()

	return stdout, stderr


def extract_info(stringblock, begin_keyword, stop_keyword, skip):
	'''
	input:
	stringblock: The file, or screen output; e.g. stdout.
	begin_keyword: The keyword user want to find, starting from this line (inclusive), the information will be extracted. If multiple begin_keywords were found, those blocks would probably all be extracted. TODO: seperate those blocks
	end_keyword: The line of ending keyword will not be extracted.
	skip: skip n line counting from begin_keyword, the information will be extracted.
	'''
	infolist = []
	appendflag = False
	startappend = False
	count_skip = 0
	for line in stringblock.splitlines():
		if line.find(begin_keyword) != -1:
			appendflag = True
			if count_skip >= skip:
				appendflag = True
		if line.find(stop_keyword) != -1:
			appendflag = False
			startappend = False
		if appendflag:
			if count_skip < skip:
				count_skip += 1
			else:
				startappend = True
		if startappend:
			infolist.append(line)

	return infolist


def custom_split(s, column_widths):
	'''
	input: 
	s: string
	column_widths: a list of column widths. This function will seperate the string accroding to the list of widths. 
	'''
	start = 0
	columns = []
	for width in column_widths:
		end = start + width
		columns.append(s[start:end])
		start = end
	return columns


class pattern_blocks:
	def __init__(self, stringblock):
		self.info_dictionary = {}
		self.endinfo_dictionary = {}
		self.allinfo_dictionary = {}
		self.stringblock = stringblock
		
		self.pattern_blocks()
	def pattern_blocks(self):
		'''
		atomatically match the pattern in the information and extract them.
		pattern should be: seperate column and find the length of string in same column different lines are more or less the same
		input:
		stringblock: The file, or screen output; e.g. stdout.
		output:
		info_dictionary:
		block title -> block content, endinfo, allinfo
		'''

		infolist = []

		count_and_compare1 = []
		count_and_compare2 = []
		pattern_flag = False
		individual_col_pattern = True
		title_line = -1
		allContent = []

		
		#for i, line in enumerate(stringblock.splitlines()):
		for line in self.stringblock:
			allContent.append(line.strip())

		for i, line in enumerate(allContent):
			tmp = line
			tmp = re.sub(r'\s+', ' ', tmp)	# 将多个连续空格替换为一个空格
			tmp = tmp.strip().split(' ')
			
			for word in tmp:
				count_and_compare2.append(len(word))
			if len(count_and_compare1) == len(count_and_compare2):
				for j, num in enumerate(count_and_compare1):
					if abs(count_and_compare1[j] - count_and_compare2[j]) > 2:
						#print(tmp)
						individual_col_pattern = False
						break


				if individual_col_pattern:
					#print("pattern found!")
					pattern_flag = True
				else:
					pattern_flag = False
			else:
				individual_col_pattern = False

			if pattern_flag:
				infolist.append(allContent[i - 1])
				
				if title_line == -1:
					title_line = i - 2


			if not individual_col_pattern:
				pattern_flag = False
				if title_line > -1 and len(infolist) > 0:
					key = allContent[title_line]
					content = infolist
					endinfo = allContent[i]
					self.info_dictionary[key] = content
					self.endinfo_dictionary[key] = endinfo
					if key not in self.allinfo_dictionary:
						self.allinfo_dictionary[key] = []
						self.allinfo_dictionary[key].append(content)
					else:
						self.allinfo_dictionary[key].append(content)
					title_line = -1
					infolist = []
					


			count_and_compare1 = copy.deepcopy(count_and_compare2)
			count_and_compare2 = []
			individual_col_pattern = True
		return self.info_dictionary, self.endinfo_dictionary, self.allinfo_dictionary

	def print_last_info(self):
		for key in self.info_dictionary.keys():
			print("title:", key)
			for line in self.info_dictionary[key]:
				print(line)

	def query_last(self, title):
		'''
		输出对应关键词title 的text block
		input:
		title: The keyword
		return:
		A list of text that contain the original information:
		'''
		if title not in self.info_dictionary:
			print("Corresponding title not found")
			return None
		else:
			return self.info_dictionary[title]

class Gaussianinfo(pattern_blocks):
	def __init__(self, stringblock):
		super().__init__(stringblock)
		self.last_xyz = []
		self.last_AtomicNum = []
		self.last_NumAtoms = 0

	def get_last_info(self):
		info = super().query_last("---------------------------------------------------------------------")
		for line in info:
			tmp = re.sub(r'\s+', ' ', line)
			tmp = tmp.strip().split(' ')
			self.last_AtomicNum.append( int(tmp[1]) )
			self.last_xyz.append( float(tmp[3]) )
			self.last_xyz.append( float(tmp[4]) )
			self.last_xyz.append( float(tmp[5]) )

		self.last_NumAtoms = len(self.last_AtomicNum)


	def print_last_xyz(self, filename=None):
		if filename is None:
			print(self.last_NumAtoms, "\n")
			for index in range(len(self.last_AtomicNum)):
				print(gc.AtomicNum_to_Symbol[self.last_AtomicNum[index]],
		  self.last_xyz[3 * index], self.last_xyz[3 * index + 1], self.last_xyz[3 * index + 2]
		  )
		#else:
