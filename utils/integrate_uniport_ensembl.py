from infra import *
import constants


def binarySearch(alist, item):
	    first = 0
	    last = len(alist)-1
	    found = False

	    while first<=last and not found:
	        midpoint = (first + last)//2
	        if alist[midpoint] == item:
	            found = True
	        else:
	            if item < alist[midpoint]:
	                last = midpoint-1
	            else:
	                first = midpoint+1

	    return found, midpoint



u2g = np.array(load_phenotype_data("uniport_to_esng.txt", phenotype_list_path=constants.DICTIONARIES_DIR))
u2g_header = u2g[0]
u2g_content = u2g[1:][u2g[1:,0].argsort(),:]


u2p = np.array(load_phenotype_data("uniport_to_esnp.txt", phenotype_list_path=constants.DICTIONARIES_DIR))
u2p_header = u2p[0]
u2p_content = u2p[1:][u2p[1:,0].argsort(),:]

integ = []
integ.append("\t".join(["Uniport Entry", "ENSP", "ENSG", "GO Terms"]))
for i, cur_entry in enumerate(u2p_content):
    print i
    row = cur_entry[:-1]

    additional = []
    found, midpoint = binarySearch(u2g_content[:,0], cur_entry[0])
    if found:
        additional = list(u2g_content[midpoint][1:])
    else:
        continue
    # for cur_entry in u2g_content:
    #     if cur_entry[0] == cur_gdc[0]:
    #         additional = list(cur_entry[1:])
    #         break
    additional = additional + [" " for x in range((len(u2g[0][1:]) - len(additional)))]
    row = list(row[0:2]) + list(additional) + list(row[2:])
    integ.append("\t".join(row))

with open(os.path.join(constants.OUTPUT_GLOBAL_DIR, "uniport_dictionary_ensg_only.txt"),"w+") as f:
    for cur in integ:
        f.write(cur+"\n")

