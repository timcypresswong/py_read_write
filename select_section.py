import argparse
import copy


def parser_args():
	parser = argparse.ArgumentParser(description='This is the program to seperate file content. If a file is ,approximately repeatative, this code recognize the repeatition and help extract the information section.')
	parser.add_argument('filename', help='input none binary file')
	parser.add_argument('--frame_number', help="input the number of frame", type = int)
	parser.add_argument('--check', help="check whether the automated recognition is correct", type = bool)
	parser.add_argument("--output", help='whether output the result in the output filename', type = str)
	args = parser.parse_args()
	return args


def run():
	args = parser_args()
	file_path = args.filename
	frame_number = args.frame_number
	check = args.check
	output = args.output
	if check is None:
		check = False
	if frame_number is None:
		frame_number = -1

	recognize_and_seperate(file_path, frame_number, check, output)
	
def recognize_and_seperate(file_path, frame_number, check, output):
	input_file = open(file_path, 'r')
	content_list = []
	charcount_list = []
	for line in input_file:
		content = line.strip()
		content_list.append(content)
		charcount_list.append(len(content))
	determined_seglength = determine_seglength(charcount_list)
	# print(determined_seglength)
	re_organized_list = [content_list[i: i + determined_seglength] for i in range(0, len(content_list), determined_seglength)]
	frame = re_organized_list[frame_number]
	if output is None:		
		for info in frame:
			print(info)
	else:
		for info in frame:
			print(info, file = output)

def determine_seglength(charcount_list):

	seglength_plus_cost = []
	for seglength in range(1, int(len(charcount_list) / 2)):
		seglength_plus_cost.append(seglength + evaluate_cost_of_list_length_every_seglength(charcount_list, seglength))
		# print(seglength, evaluate_cost_of_list_length_every_seglength(charcount_list, seglength))
	min_cost = min(seglength_plus_cost)
	determined_seglength = seglength_plus_cost.index(min_cost) + 1
	return determined_seglength

def evaluate_discrete_function_difference(intlist1, intlist2) -> int:
    maxlength = max( len(intlist1), len(intlist2))
    smalllist = None
    biglist = None
    if len(intlist1) < maxlength:
        smalllist = copy.deepcopy(intlist1)
        biglist = copy.deepcopy(intlist2)
    else:
        smalllist = copy.deepcopy(intlist2)
        biglist = copy.deepcopy(intlist1)
    for i in range(len(smalllist), maxlength):
        smalllist.append(0)

    sum_diff = 0
    for i in range(maxlength):
        diff = abs(smalllist[i] - biglist[i])**2
        sum_diff += diff
    return sum_diff


def evaluate_cost_of_list_length_every_seglength(intlist: list, seglength: int) -> int:
    '''
    Docstring for evaluate_cost_of_list_length_every_seglength

    :param intlist: the entire list of length of character
    :type intlist: list
    :param seglength: number of line to consider a segtion
    :type seglength: int
    :return: the value of cost function
    :rtype: int
    '''
    re_organized_list = [intlist[i: i + seglength] for i in range(0, len(intlist), seglength)]
    num_fraction = len(re_organized_list)
    cost = 0
    for i in range(num_fraction):
        for j in range(i+1, num_fraction):
            cost += evaluate_discrete_function_difference(re_organized_list[i], re_organized_list[j])
    return cost /( num_fraction - 1)


if __name__ == '__main__':
	run()

