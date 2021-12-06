import os
import math


def operation_on_cells(file1, file2, operation, output_file='output.tsv',\
	sep='\t', eol='\n', col_header=True, row_header=False, decimals=3):
	"""file1 and file2 being two files of same dimensions, applies the 
	chosen operation (+-*/)	on every cell of file1 with its homolog in 
	file2."""

	assert operation in '+-*/'
	assert output_file != ''

	max_result = 0

	error = False
	output_content = ''
	with open(file1, 'r') as f1:
		with open(file2, 'r') as f2:
			if col_header:
				# Add column headers (should be identical in both files)
				output_content += f1.readline()
				f2.readline()

			# Calculate results for all lines
			line1 = f1.readline()
			line2 = f2.readline()
			while (len(line1), len(line2)) != (0, 0) and not error:
				line1_list = line1.rstrip(eol).split(sep)
				line2_list = line2.rstrip(eol).split(sep)

				if row_header:
					line_results = [line1_list[0]]
				else:
					line_results = []

				for i in range(1, len(line1_list)-1):
					try:
						line1_list[i] = float(line1_list[i])
						line2_list[i] = float(line2_list[i])
					except ValueError:
						print('Some values were not numbers')
						error = True
						break

					if operation == '+':
						result = line1_list[i]+line2_list[i]
					elif operation == '-':
						result = line1_list[i]-line2_list[i]
					elif operation == '*':
						result = line1_list[i]*line2_list[i]
					elif operation == '/':
						result = line1_list[i]/line2_list[i]
					else:
						print('Undefined operation')
						error = True
						break

					if result > max_result:
						max_result = result

					format_result = '{:.'+str(decimals)+'f}'
					result = format_result.format(result)
					line_results.append(result)
				
				output_content += sep.join(line_results)+eol

				# average = 0
				# for result in line_results[1:]:
				# 	average += float(result)
				# average /= len(line_results[1:])
				# print(line_results[0], average)

				line1 = f1.readline()
				line2 = f2.readline()

	print('max_result:', max_result)

	with open(output_file, 'w') as f:
		f.write(output_content)
		print(output_file, 'written')

	return True


if __name__ == '__main__':
	# mydir = os.path.dirname(__file__)
	# dN_file = mydir+'/'+'example2_dN.count'
	# dS_file = mydir+'/'+'example2_dS.count'

	operation_on_cells('example2_dN.count', 'example2_dS.count', \
		operation='/', output_file='example2_omegas.count', \
		col_header=True, row_header=True, decimals=4)

	# input()
