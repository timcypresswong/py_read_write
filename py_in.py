import subprocess
import os
import pandas as pd


def find_path():
	script_path = os.path.abspath(__file__)
	script_dir = os.path.dirname(script_path)
	pathrc = script_dir + "/.appspathrc"
	pathdf = pd.read_csv(pathrc, sep='\t', header=None)
	pathdict = pd.Series(pathdf[1].values, index=pathdf[0]).to_dict()
	return pathdict


def call_multiwfn(file_path, commands):
#	multiwfn_path = "/home/dunbo/software/compchem/Multiwfn/Multiwfn_3.8_dev_bin_Linux_noGUI/Multiwfn_noGUI"
	multiwfn_path = find_path()['multiwfn']

	fchk = open(file_path, 'rb')

	process = subprocess.Popen([multiwfn_path, file_path], stdin=subprocess.PIPE,
	stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

	stdout, stderr = process.communicate(input=commands)
	fchk.close()

	return stdout, stderr


def extract_info(stringblock, begin_keyword, stop_keyword, skip):
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
	start = 0
	columns = []
	for width in column_widths:
		end = start + width
		columns.append(s[start:end])
		start = end
	return columns
