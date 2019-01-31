from infra import *
import constants




go_dict = {}

def load_dict():
	u_d = np.array(load_phenotype_data("go_dictionary.txt", phenotype_list_path=constants.DICTIONARIES_DIR))
	u_d_content = u_d[1:][u_d[1:, 0].argsort(), :]

	for i, cur_entry in enumerate(u_d_content):
		if go_dict.has_key(cur_entry[0]):
			go_dict[cur_entry[0]]["Uniport Entry"].add(cur_entry[1])
			go_dict[cur_entry[0]]["ENSP"].add(cur_entry[2])
			go_dict[cur_entry[0]]["ENSG"].add(cur_entry[3])
		else:
			go_dict[cur_entry[0]] = {"Uniport Entry" : set(cur_entry[1].split("; ")), "ENSP" : set(cur_entry[2].split("; ")), "ENSG" : set(cur_entry[3].split("; "))}
	return go_dict

def get_go_dict():
	ensg_dict = load_dict()
	return  ensg_dict

