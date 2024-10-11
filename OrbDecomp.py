import py_in
import argparse

def parser_args():
	parser = argparse.ArgumentParser(description='This is the program for assigning MO to AO.')
	parser.add_argument('filename', help='input formchk file')
	parser.add_argument('--orblist', help='input a list of orbital', nargs='+', type=str)
	args = parser.parse_args()
	return args




def run():
	args = parser_args()
	file_path = args.filename
	orbital_list = args.orblist

	MullikenMainOrb(file_path, orbital_list)


def MullikenMainOrb(file_path, orbital_list):

	atomlist = []
	typelist = []
	locallist = []

	for orb in orbital_list:
		commands = f"8\n1\n{orb}\n0\n-10\nq\n"
		stdout, stderr = py_in.call_multiwfn(file_path, commands)

		info = py_in.extract_info(stdout, "Basis Type    Atom    Shell      Local       Cross term        Total",  
									"Sum up those listed above:", 1)
		for line in info:
			columns = py_in.custom_split(line, [7, 5, 11, 7, 13, 13, 13])
			atomlist.append(columns[2])
			typelist.append(columns[1])
			locallist.append( float(columns[4].split()[0]) / 100 )

		max_value = max(locallist)
		max_position = locallist.index(max_value)

		print("orbital label:", orb)
		print("main composition:", atomlist[max_position], typelist[max_position] , "with value:", max_value )
		print("summary:")
		for i,x in enumerate(atomlist):
			print(typelist[i], atomlist[i], locallist[i])


		atomlist = []
		typelist = []
		locallist = []


#file_path = "/home/dunbo/software/compchem/Multiwfn/script/graphene.fch"
#orbital_list = ['h', 'h-1', 'h-2']

if __name__ == '__main__':
	run()
