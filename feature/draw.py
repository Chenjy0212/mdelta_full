import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--json', '-j', help='[path/filename]Result file after running mdelta(in JSON format)', required=True)
parser.add_argument('--output', '-o', help='[path] Output path, eg:\'/home/username\'')
parser.add_argument('-nt','--Name2TypeFile',type=str, nargs='?', default='', help='[path/filename] List of correspondance between tip name and cell type for cell lineage tree #1.')
parser.add_argument('-nt2','--Name2TypeFile2',type=str, nargs='?', default='', help='[path/filename] List of correspondance between tip name and cell type for cell lineage tree #2.')
parser.add_argument('-xsd','--XScoreDictFile',type=str, nargs='?', default='', help='[path/filename] A comma-delimited text file used to determine similarity scores between cells. If there are exactly three columns, they will be interpreted as (1) the cell (name or type) in Tree #1, (2) the cell in Tree #2, and (3) the similarity score. If otherwise, the first column will be interpreted as the cell (name or type) and the remaining columns as features of the cell (e.g. expression of a gene). The similarity scores will be estimated between all pairs of cells based on the Euclidean distance calculated using all the features. Overrides `-ma` and `-mi`.')
parser.add_argument('-ma','--mav',type=float, default=2., help=' [float] Default=2.')
# 0$result_file 1$output 2$name2type1 3$name2type2 4$scoredict 5$mav

args = parser.parse_args()
print(args)

# os.system('Rscript /home/ee_while/JOB230605/match_tree.r {0[0]} {0[4]} {0[1]} {0[5]}'.format(args.json))
# print('match_tree ok!!!')
    
# os.system('python /home/ee_while/JOB230605/network.py {0[0]} {0[1]}'.format(args.list))
# print('network ok!!!')

# os.system('Rscript /home/ee_while/JOB230605/densitree.r {0[0]} {0[2]} {0[3]} {0[1]}'.format(args.list))
# # os.system('Rscript /home/ee_while/JOB230605/densitree.r {0[0]} {0[1]}'.format(args.list))
# print('densitree ok!!!')