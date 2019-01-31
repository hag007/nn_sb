from infra import *
import constants




ensg_dict = None

def get_ensg_dict():
	global ensg_dict
	if ensg_dict is None:
		ensg_dict = load_dict()
	return ensg_dict

def load_dict():
	global ensg_dict
	ensg_dict = {}
	u_d = np.array(load_phenotype_data("uniport_dictionary_ensg_only.txt", phenotype_list_path=constants.DICTIONARIES_DIR))
	u_d_header = u_d[0]
	u_d_content = u_d[1:][u_d[1:, 0].argsort(), :]

	for i, cur_entry in enumerate(u_d_content):
		if ensg_dict.has_key(cur_entry[0]):
			# ensg_dict[cur_entry[0]]["ENSP"].add(cur_entry[1])
			# ensg_dict[cur_entry[0]]["ENSG"].add(cur_entry[2])
			ensg_dict[cur_entry[0]]["GO Terms"].update(cur_entry[3].split("; "))
		else:
			ensg_dict[cur_entry[0]] = {"ENSP" : set([cur_entry[1]]), "ENSG" : set([cur_entry[2]]), "GO Terms" : set(cur_entry[3].split("; "))}
	return ensg_dict



