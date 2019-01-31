import simplejson as json
from utils.ensembl2gene_symbol import g2e_convertor

js=json.loads(file("/home/hag007/Downloads/Parsimonious Composite Network (PCNet).cx").read())


edges=js[7]
nodes=js[8]

dict_nodes={}
for cur in nodes["nodes"]:
    id=cur["@id"]
    n=cur["n"]
    dict_nodes[id]=n


output="ID_interactor_A\tppi\tID_interactor_B\n"
for cur in edges["edges"]:
    s=g2e_convertor(dict_nodes[cur["s"]])[0] if len(g2e_convertor(dict_nodes[cur["s"]])) > 0 else ""
    t=g2e_convertor(dict_nodes[cur["t"]])[0] if len(g2e_convertor(dict_nodes[cur["t"]])) > 0 else ""
    if s!="" and "R" not in s and t!="" and "R" not in t:
        output+="{}\tppi\t{}\n".format(s ,t)


file("/media/hag007/Data/bnet/networks/PCNet.sif", 'w+').write(output[:-2])


