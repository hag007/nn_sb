import sys
sys.path.insert(0, '../')
from infra import *
import constants





u_d = np.array(load_phenotype_data("esnp_dictionary.txt", phenotype_list_path=constants.DICTIONARIES_DIR))
u_d_header = u_d[0]
u_d_content = u_d[1:][u_d[1:,0].argsort(),:]

ensg_dict = {}
integ = []
integ.append("\t".join(["GO Term", "ENSP", "Uniport Entry", "ENSG"]))
for i, cur_entry in enumerate(u_d_content):
	go_terms = cur_entry[3].split("; ")
	for cur_go_term in go_terms:
		cur_go_term = cur_go_term.strip()
		print i
		if ensg_dict.has_key(cur_go_term):
			ensg_dict[cur_go_term]["ENSP"].add(cur_entry[0])
			ensg_dict[cur_go_term]["Uniport Entry"].add(cur_entry[1])
			ensg_dict[cur_go_term]["ENSG"].add(cur_entry[2])
		else:
			ensg_dict[cur_go_term] = {"ENSP" : set([cur_entry[0]]), "Uniport Entry" : set([cur_entry[1]]), "ENSG" : set([cur_entry[2]])}

for k, v in ensg_dict.iteritems():
	row = "{}\t{}\t{}\t{}".format(
		k,
		"; ".join(v["ENSP"]),
		"; ".join(v["Uniport Entry"]),
		"; ".join(v["ENSG"]))
	integ.append(row)


with open(os.path.join(constants.OUTPUT_GLOBAL_DIR, "go_dictionary.txt"),"w+") as f:
    for cur in integ:
        f.write(cur+"\n")

