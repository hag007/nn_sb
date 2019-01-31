import constants
import os 

go_dict = None


def get_go_names(GO_ids):
    GO_names = []
    global go_dict
    if go_dict is None:
        go_dict = {}
        f = file(os.path.join(constants.GO_DIR, 'go-basic.obo'))
        parsed = f.read().split("\n\n")
        for cur_obo in parsed[1:]:
            if cur_obo.split("\n")[0] != "[Term]": continue
            go_dict[cur_obo.split("\n")[1][4:]] = cur_obo.split("\n")[2][6:]

    for cur_id in GO_ids:
        if cur_id in go_dict:
            GO_names.append(go_dict[cur_id])
        else:
            GO_names.append(cur_id)

    # print "\n".join(GO_names)
    return GO_names

