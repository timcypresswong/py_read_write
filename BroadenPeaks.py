import py_in
import numpy as np
import argparse
import pandas as pd
import os

def parser_args():
	parser = argparse.ArgumentParser(description='This is the program that read a controlling yaml to broaden columnic data.')
	parser.add_argument('filename', help='input yaml file')

	args = parser.parse_args()
	return args




def run():
	args = parser_args()
	yaml_path = args.filename
	inputyaml = py_in.yaml_base(yaml_path)
	dir_path = None
	outputfileflag = False
	abs_here = os.getcwd()
	abs_yaml = os.path.dirname(args.filename)
	default_dir_paths = [abs_here, abs_yaml]												# search the input file either in "." or in the 
	#inputyaml.print_yaml_info()															# same path of yaml file
	if not inputyaml.exist_key_1st("filename"):
		raise KeyError("key filename not in the yaml input")
	if not inputyaml.exist_key_1st("kernel"):
		raise KeyError("Please choose a kernel function")
	

	Peakfilename = inputyaml.data["filename"]
	kernel = inputyaml.data["kernel"]

	if inputyaml.exist_key_1st("dir_path"):
		dir_path = inputyaml.data["dir_path"]
		if isinstance(Peakfilename, str):
			tmp_file = os.path.join(dir_path, Peakfilename)
			if os.path.exists(tmp_file):
				Peakfilename = tmp_file
			else:
				raise FileNotFoundError(f"file {tmp_file} does not exists")
		if isinstance(Peakfilename, list):
			absPeakfilename = []
			for name in Peakfilename:
				tmp_file = os.path.join(dir_path, name)
				if os.path.exists(tmp_file):					
					absPeakfilename.append( os.path.join(dir_path, tmp_file) )
				else:
					raise FileNotFoundError(f"file {tmp_file} does not exists")
			Peakfilename = absPeakfilename

	else:
		found_flag = False
		
		for _path in default_dir_paths:
			if isinstance(Peakfilename, str):
				tmp_file = os.path.join(_path, Peakfilename)
				if os.path.exists(tmp_file):
					Peakfilename = tmp_file
					found_flag = True
					break
			if isinstance(Peakfilename, list):
				absPeakfilename = []
				for name in Peakfilename:
					tmp_file = os.path.join(_path, name)
					if os.path.exists(tmp_file):					
						absPeakfilename.append( tmp_file )
						found_flag = True
						dir_path = _path
					else:
						found_flag = False
						absPeakfilename = []
						dir_path = _path
						break
				if absPeakfilename != []:
					Peakfilename = absPeakfilename

		if found_flag == False:
			raise FileNotFoundError(f"file {Peakfilename} does not exists in the default paths: {default_dir_paths}")

	
	
	resolution = 0.01
	if inputyaml.exist_key_1st("resolution") and inputyaml.data["resolution"] is not None:
		resolution = float(inputyaml.data["resolution"])
	

	x = 1																					# broadening x axis
	if inputyaml.exist_key_1st("x") and inputyaml.data["x"] is not None:
		x = int(inputyaml.data["x"])

	y = 2
	if inputyaml.exist_key_1st("y") and inputyaml.data["y"] is not None:
		y = int(inputyaml.data["y"])


	if inputyaml.exist_key_1st("output") and inputyaml.data["output"] == True:
		outputfileflag = True

	'''
	parameters for Lorentz
	'''
	FWHM = 0.05
	if kernel == "Lorentz":
		if inputyaml.exist_key_1st("FWHM") and inputyaml.data["FWHM"] is not None:
			FWHM = float(inputyaml.data["FWHM"])

		xb = None
		yb = None
		if isinstance(Peakfilename, list):
			for filename in Peakfilename:
				xb,yb = Broaden_Lorentz(filename, x,y,resolution, FWHM)
				printtwocolumns(xb,yb,outputfileflag, newname(filename,"_Lorentz"))
		if isinstance(Peakfilename, str):
			xb, yb = Broaden_Lorentz(Peakfilename, x,y,resolution, FWHM)
			printtwocolumns(xb,yb,outputfileflag, newname(Peakfilename,"_Lorentz"))
		
def printtwocolumns(x, y, outputfileflag, path):
	if outputfileflag == False:
		for i in range(len(x)):
			print(f"{x[i]}\t{y[i]}")
		print("")
	else:
		file = open(path, 'w')
		for i in range(len(x)):
			print(f"{x[i]}\t{y[i]}", file = file)


def newname(path, newidentifier ):
	body = path.rsplit('.')[0]
	suffix = path.rsplit('.')[1]
	return body + newidentifier + '.' + suffix



def Broaden_Lorentz(filename, x, y, resolution, FWHM):
	'''
	input:
	filename: the filepath
	x: the column idx for data on x axis (1-based)
	y: the column idx for data on y axis (1-based)
	resolution: the interval between points of the curve
	FWHM: the broadening factor

	return:
	xrange: broadened data on x axis
	y: braodened data on y axis
	'''
	Peaks = pd.read_csv(filename, sep='\t', comment='#', header=None)
	colx = np.array(Peaks[x-1].tolist())
	coly = np.array(Peaks[y-1].tolist())
	xPeakmin = np.min(colx)
	xPeakmax = np.max(colx)
	xmin = xPeakmin - 50 * FWHM
	xmax = xPeakmax + 50 * FWHM
	xrange = np.array(np.linspace(xmin,xmax, int((xmax - xmin) / resolution) ))
	y = LorentzFuncMultiPeaks(xrange, colx, coly, FWHM)
	return xrange, y

def LorentzFunc(x,x0,FWHM):
	y = np.zeros(len(x))
	for i , xi in enumerate(x):
		#print("x0:", x0, "xi:", xi)
		y[i] = 1 / np.pi * (FWHM / 2) / ((xi-x0)**2 + (FWHM / 2)**2)
	return y

def LorentzFuncMultiPeaks(x,xs,ys,FWHM):
	y = np.zeros(len(x))
	for i, peakx in enumerate(xs):
		y += ys[i] * LorentzFunc(x, peakx, FWHM)
	return y



if __name__ == '__main__':
	run()
