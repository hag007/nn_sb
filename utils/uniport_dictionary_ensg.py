from infra import *
import constants





u_d = np.array(load_phenotype_data("esnp_dictionary.txt", phenotype_list_path=constants.DICTIONARIES_DIR))
u_d_header = u_d[0]
u_d_content = u_d[1:][u_d[1:,0].argsort(),:]

ensg_dict = {}
integ = []
integ.append("\t".join(["ENSG", "ENSP", "Uniport Entry", "GO Terms"]))
for i, cur_entry in enumerate(u_d_content):
	print i
	if ensg_dict.has_key(cur_entry[2]):
		ensg_dict[cur_entry[2]]["ENSP"].add(cur_entry[0])
		ensg_dict[cur_entry[2]]["Uniport Entry"].add(cur_entry[1])
		ensg_dict[cur_entry[2]]["GO Terms"].update(cur_entry[3].split("; "))
	else:
		ensg_dict[cur_entry[2]] = {"ENSP" : set([cur_entry[0]]), "Uniport Entry" : set([cur_entry[1]]), "GO Terms" : set(cur_entry[3].split("; "))}

for k, v in ensg_dict.iteritems():
	row = "{}\t{}\t{}\t{}".format(
		k,
		"; ".join(v["ENSP"]),
		"; ".join(v["Uniport Entry"]),
		"; ".join(v["GO Terms"]))
	integ.append(row)


with open(os.path.join(constants.OUTPUT_GLOBAL_DIR, "esng_dictionary.txt"),"w+") as f:
    for cur in integ:
        f.write(cur+"\n")

