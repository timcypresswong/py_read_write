import py_in
import argparse

def parser_args():
	parser = argparse.ArgumentParser(description='This is the program to extract the pattern information from a file.')
	parser.add_argument('filename', help='input none binary file')
	parser.add_argument('--frame_number', help="input the number of frame", type = int)
	parser.add_argument("--all_frame", help='whether or not output all frame', type = bool)
	parser.add_argument("--output", help='whether output the result in the output filename', type = str)
	#TODO: also able to handle stdout
	args = parser.parse_args()
	return args




def run():
	args = parser_args()
	file_path = args.filename
	frame_number = args.frame_number
	all_frame = args.all_frame
	output = args.output
	if all_frame is None:
		all_frame = False
	if frame_number is None:
		frame_number = -1

	printGaussianlog(file_path, frame_number, all_frame, output)



def printGaussianlog(file_path, frame_number, all_frame, output):
	input_file = open(file_path, 'r')
	xyzinfo = py_in.Gaussianinfo(input_file)
	xyzinfo.get_all_info()

	# xyzinfo.print_all_xyz()
	if all_frame:
		xyzinfo.print_all_xyz(output)
		return 0
	else:
		xyzinfo.print_frame_xyz(frame_number, output)

'''	input_file.close()
	input_file = open(file_path, 'r')
	cleverinfo = py_in.pattern_blocks(input_file)
	infoblock = cleverinfo.query_all_by_key("---------------------------------------------------------------------")'''
	#print(infoblock)



if __name__ == '__main__':
	run()
