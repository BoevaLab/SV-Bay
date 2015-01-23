f = open('/Users/dasha/PhD/results_3_tools/links_probabilities.txt','r')
links_models=dict()
for line in f:
	if 'Best' in line:
		line = line.split()
		links_models[line[4]] = dict()
		links_models[line[4]]['A1'] = line[11][1:-1]
		links_models[line[4]]['B2'] = line[14][:-1]
f.close()
valid = open('/Users/dasha/PhD/results_3_tools/valid_links/all_valid_pe.txt','r')
valid_new = open('/Users/dasha/PhD/results_3_tools/all_valid_wt_model_pe.txt','w')
for line in valid:
	print line
	line_sp = line.split(';')
	new_line=line[:-1]+';'+links_models[line_sp[0]]['A1']+';'+links_models[line_sp[0]]['B2']+'\n'
	print new_line
	valid_new.write(new_line)
valid.close()
valid_new.close()

