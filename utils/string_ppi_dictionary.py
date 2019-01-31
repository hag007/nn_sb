from infra import *
import constants




string_ppi_dict = {}

def load_dict():
	u_d = np.array(load_phenotype_data("string_ppi_total.txt", phenotype_list_path=constants.DICTIONARIES_DIR))
	u_d_content = u_d[1:]

	for i, cur_entry in enumerate(u_d_content):
		edge = "{}={}".format(cur_entry[0][5:], cur_entry[1][5:])
		string_ppi_dict[edge] = cur_entry[2]
	return string_ppi_dict

def get_string_ppi_dict():
	string_ppi_dict = load_dict()
	return string_ppi_dict

