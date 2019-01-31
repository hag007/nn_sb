from infra import *
import constants





u_d = np.array(load_phenotype_data("uniport_dictionary_ensg_only.txt", phenotype_list_path=constants.DICTIONARIES_DIR))
u_d_header = u_d[0]
u_d_content = u_d[1:][u_d[1:,0].argsort(),:]

integ = []
integ.append("\t".join(["ENSP", "Uniport Entry", "ENSG", "GO Terms"]))
for i, cur_entry in enumerate(u_d_content):
	print i
	ensp_ids = cur_entry[1].split(",")
	for cur in ensp_ids:
		integ.append("\t".join([cur] + [cur_entry[0]] + list(cur_entry[2:])))


with open(os.path.join(constants.OUTPUT_GLOBAL_DIR, "esnp_dictionary.txt"),"w+") as f:
    for cur in integ:
        f.write(cur+"\n")

