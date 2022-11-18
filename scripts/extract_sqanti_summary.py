import argparse
import re

parser = argparse.ArgumentParser(description='Extract textual summary to stdout from a SQANTI report')
parser.add_argument('html', metavar='HTML',
                    help='SQANTI report HTML file')
args = parser.parse_args()

html_path = args.html


data_prog = re.compile(r".*\"data\":\[(?P<data>.*)\],\"container\".*")
list_prog = re.compile(r"\[(?P<list>[^\[\]]*)\]")

def extract_table_columns(line):
    data_result = data_prog.match(line)
    data_str = data_result.group('data')
    list_result = list_prog.findall(data_str)
    
    columns = list()
    for list_str in list_result[1:]:
        columns.append([i.strip('"').replace('\\n', ' ') for i in list_str.split(',')])
    
    return columns

def print_columns(prefix, columns):
    for line in zip(*columns):
        # only print the first two columns
        print(prefix, line[0], line[1], sep='\t')

seen_gene = False
seen_isoform = False
seen_sj = False        

print('classification', 'category', 'count', sep='\t')

with open(html_path) as fr:
    for line in fr:
        if seen_gene and seen_isoform and seen_sj:
            break
        else:
            line = line.strip()
            if not seen_gene and line.casefold() == "<h3>Gene classification</h3>".casefold():
                fr.readline() # heading line
                cols = extract_table_columns(fr.readline())
                print_columns('gene', cols)
                seen_gene = True
            elif not seen_isoform and line.casefold() == "<h3>Isoform Classification</h3>".casefold():
                fr.readline() # heading line
                cols = extract_table_columns(fr.readline())
                print_columns('isoform', cols)
                seen_isoform = True
            elif not seen_sj and line.casefold() == "<h3>Splice Junction Classification</h3>".casefold():
                fr.readline() # heading line
                cols = extract_table_columns(fr.readline())
                print_columns('junction', cols)
                seen_sj = True

