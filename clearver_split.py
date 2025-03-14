import py_in
import argparse

def parser_args():
	parser = argparse.ArgumentParser(description='This is the program to extract the pattern information from a file.')
	parser.add_argument('filename', help='input none binary file') #TODO: also able to handle stdout
	args = parser.parse_args()
	return args




def run():
	args = parser_args()
	file_path = args.filename

	printtest(file_path)



def printtest(file_path):
	input_file = open(file_path, 'r')
	#cleverinfo = py_in.pattern_blocks(input_file)
	#infoblock = cleverinfo.query_last("---------------------------------------------------------------------")
	#print(infoblock)
	xyzinfo = py_in.Gaussianinfo(input_file)
	xyzinfo.get_last_info()
	xyzinfo.print_last_xyz()



if __name__ == '__main__':
	run()
