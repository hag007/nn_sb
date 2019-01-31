###
# Download the GO database and gene2go annotations
# Write out file with the GO hierarchy for human:
# <child> \t <parent> \t <child-parent link>
# child to parent link can be: gene, is_a, part_of, [positively_/negatively_]regulates
###
import time
from goatools import obo_parser
from  goatools.associations import read_ncbi_gene2go
from goatools.base import download_ncbi_associations
import constants

import wget, os

def fetch_go_hierarcy():

    obo_file_location = os.path.join(constants.GO_DIR,constants.GO_FILE_NAME)
    if not os.path.exists(os.path.join(constants.GO_DIR,constants.GO_FILE_NAME)):
            wget.download(constants.GO_OBO_URL, os.path.join(constants.GO_DIR,constants.GO_FILE_NAME))

    go = obo_parser.GODag(obo_file_location,optional_attrs=['relationship']) # also use

    print "Downloading gene-GO associations"
    association_file_location = os.path.join(constants.GO_DIR,constants.GO_ASSOCIATION_FILE_NAME)
    if not os.path.exists(association_file_location):
            wget.download(constants.GO_ASSOCIATION_GENE2GEO_URL, os.path.join(constants.GO_DIR,constants.GO_ASSOCIATION_FILE_NAME))

    print "Loading gene-GO associations"
    # gene2go = download_ncbi_associations(obo_file_location) - why does this line needed?
    go2geneids_human = read_ncbi_gene2go(association_file_location, taxids=[9606], go2geneids=True)


    print "Writing out GO child-parent links"
    if not os.path.exists(constants.OUTPUT_GLOBAL_DIR):
            os.makedirs(constants.OUTPUT_GLOBAL_DIR)

    out_fname = "go_output_{}_{}.txt".format(constants.CANCER_TYPE, time.time())
    genes = []
    isa = []
    relship = []
    with open(os.path.join(constants.OUTPUT_GLOBAL_DIR,out_fname),'w') as o:
        for goid in go2geneids_human.keys():
            if not go.has_key(goid):
                print "GO obo file does not contain {}".format(goid)
                continue
            entry = go[goid]
            for gene in go2geneids_human[entry.id]:
                genes.append((str(gene), entry.id))
                o.write("{}\t{}\t{}\n".format("genes", *genes[-1]))
            children = entry.children
            for c in children:
                isa.append((c.id, entry.id))
                o.write("{}\t{}\t{}\n".format("is a", *isa[-1]))
            rels = entry.relationship_rev
            for rtype in rels.keys():
                rs = rels[rtype]
                for r in rs:
                    relship.append((rtype, r.id, entry.id))
                    o.write("{}\t{}\t{}\n".format(rtype, *relship[-1]))

    return (genes, isa, relship)




fetch_go_hierarcy()