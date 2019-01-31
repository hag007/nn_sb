from infra import *
import constants




prize_dict = {}

def load_dict():
	u_d = np.array(load_phenotype_data("prize_dictionary.txt", phenotype_list_path=constants.DICTIONARIES_DIR))
	u_d_header = u_d[0]
	u_d_content = u_d[1:][u_d[1:, 0].argsort(), :]

	for i, cur_entry in enumerate(u_d_content):
		prize_dict[cur_entry[0]]=cur_entry[1]
	return prize_dict

def get_prize_dict():
	prize_dict = load_dict()
	return prize_dict

