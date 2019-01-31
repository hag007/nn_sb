"""Print a GO term's lower-level hierarchy."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
sys.path.insert(0, '../')
import collections as cx
from goatools.godag.consts import Consts
from goatools.gosubdag.go_paths import GoPaths
import constants
"""driver imports"""
import os
import timeit
import datetime
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.rpt.write_hierarchy import WrHierGO
""""gene to go mapper imports"""
import time
from goatools import obo_parser
from  goatools.associations import read_ncbi_gene2go
from goatools.base import download_ncbi_associations
import constants
from utils.ensembl2entrez import get_entrez2ensembl_dictionary
import wget, os

class WrHierGO(object):
    """Write hierarchy object."""

    kws_dct = set(['max_indent'])
    kws_set = set(['no_indent', 'concise'])
    consts = Consts()

    def __init__(self, gosubdag, **kws):
        self.gosubdag = gosubdag  # GoSubDag arg, children=True, must be used
        self.usrdct = {k:v for k, v in kws.items() if k in kws}
        self.usrset = set([k for k, v in kws.items() if k in kws and v])
        # ' {NS} {dcnt:6,} L{level:02} D{depth:02} {D1:5} {GO_name}'

    def prt_hier_all(self, prt=sys.stdout):
        """Write hierarchy for all GO Terms in obo file."""
        # Print: [biological_process, molecular_function, and cellular_component]
        gos_printed = set()
        for goid in ['GO:0008150', 'GO:0003674', 'GO:0005575']:
            gos_printed.update(self.prt_hier_down(goid, prt))
        return gos_printed

    def prt_hier_down(self, goid, prt=sys.stdout):
        """Write hierarchy for all GO IDs below GO ID in arg, goid."""
        obj = _WrHierPrt(self, prt)
        obj.prt_hier_rec(goid)
        return obj.gos_printed

    def ext_hier_down(self, goid, prt=sys.stdout):
        """Write hierarchy for all GO IDs below GO ID in arg, goid."""

        obj = _WrHierPrt(self, prt)
        obj.ext_hier_rec(goid)
        # prt.write("VERTICES:\n")
        # for k, v in obj.vertices.iteritems():
        #     prt.write(str("{} - {}\n".format(k,v)))

        # prt.write("EDGES:\n")
        # for k, v in obj.edges.iteritems():
        #     prt.write(str("{} - {}\n".format(k,v)))

        return obj

    def prt_hier_up(self, goids, prt=sys.stdout):
        """Write hierarchy for all GO IDs below GO ID in arg, goid."""
        go2goterm_all = {go:self.gosubdag.go2obj[go] for go in goids}
        objp = GoPaths()
        gos_printed = set()
        for namespace, go2term_ns in self._get_namespace2go2term(go2goterm_all).items():
            go_root = self.consts.NAMESPACE2GO[namespace]
            goids_all = set()  # GO IDs from user-specfied GO to root
            for goid, goterm in go2term_ns.items():
                goids_all.add(goid)
                paths = objp.get_paths_from_to(goterm, goid_end=None, dn0_up1=True)
                goids_all.update(set(o.id for p in paths for o in p))
            # Only include GO IDs from user-specified GO to the root
            if 'include_only' not in self.usrdct:
                self.usrdct['include_only'] = set()
            self.usrdct['include_only'].update(goids_all)
            # Mark the user-specfied GO term
            if 'go_marks' not in self.usrdct:
                self.usrdct['go_marks'] = set()
            self.usrdct['go_marks'].update(go2term_ns.keys())
            obj = _WrHierPrt(self, prt)  # , goids_all, set(go2term_ns.keys()))
            gos_printed.update(obj.gos_printed)
            obj.prt_hier_rec(go_root)
        return gos_printed

    @staticmethod
    def _get_namespace2go2term(go2terms):
        """Group GO IDs by namespace."""
        namespace2go2term = cx.defaultdict(dict)
        for goid, goterm in go2terms.items():
            namespace2go2term[goterm.namespace][goid] = goterm
        return namespace2go2term


# pylint: disable=too-many-instance-attributes,too-few-public-methods
class _WrHierPrt(object):
    """Print GO hierarchy."""

    def __init__(self, obj, prt=sys.stdout):
        self.gosubdag = obj.gosubdag
        self.max_indent = obj.usrdct.get('max_indent')
        self.include_only = obj.usrdct['include_only'] if 'include_only' in obj.usrdct else None
        self.go_marks = obj.usrdct['go_marks'] if 'go_marks' in obj.usrdct else set()
        self.concise_prt = 'concise' in obj.usrset
        self.indent = 'no_indent' not in obj.usrset
        self.go2geneids = obj.usrdct.get('go2geneids')
        # vars
        self.prt = prt
        self.edges = {}
        self.vertices = {}
        self.gos_printed = set()
        self.prtfmt = self._init_prtfmt()
        self.dash_len = obj.usrdct.get('dash_len', 6) + 12

    def ext_hier_rec(self, goid, depth=1):
        """Write hierarchy for a GO Term record and all GO IDs down to the leaf level."""
        ntgo = self.gosubdag.go2nt[goid]
        ntobj = self.gosubdag.go2obj[goid]
        # Shortens hierarchy report by only printing the hierarchy
        # for the sub-set of user-specified GO terms which are connected.
        if self.include_only and goid not in self.include_only:
            return
        nrp = self.concise_prt and goid in self.gos_printed
        if self.vertices.has_key(goid):
            self.vertices[goid]["weight"] += 1
            self.vertices[goid]["depth"].append(depth)
        else:
            self.vertices[goid] = {"name": ntgo.GO_name,
                                   "weight": 0,
                                   "NS": ntgo.NS,
                                   "depth": [depth],
                                   "L" : ntgo.level,
                                   "D" : ntgo.depth,
				   "obj" : ntobj,
				   "n_children" : len(ntobj.children)}

        self.gos_printed.add(goid)
        # Do not extract hierarchy below this turn if it has already been printed
        # if nrp:
        #     return
        depth += 1
        if self.max_indent is not None and depth > self.max_indent:
            return
        for child in ntobj.children:
            if self.go2geneids.has_key(child.id) or True:
                if self.edges.has_key("{}={}".format(goid, child.id)):
                    self.edges["{}={}".format(goid, child.id)]["weight"] += 1
                else:
                    self.edges["{}={}".format(goid, child.id)] = {"weight" : 0}
                self.ext_hier_rec(child.id, depth)

    def prt_hier_rec(self, goid, depth=1):
        """Write hierarchy for a GO Term record and all GO IDs down to the leaf level."""
        ntgo = self.gosubdag.go2nt[goid]
        ntobj = self.gosubdag.go2obj[goid]
        # Shortens hierarchy report by only printing the hierarchy
        # for the sub-set of user-specified GO terms which are connected.
        if self.include_only and goid not in self.include_only:
            return
        nrp = self.concise_prt and goid in self.gos_printed
        if self.go_marks:
            self.prt.write('{} '.format('>' if goid in self.go_marks else ' '))

        # '-' is default character indicating hierarchy level
        # '=' is used to indicate a hierarchical path printed in detail previously.
        dashgo = self._str_dashgoid(ntgo, depth, not nrp or not ntobj.children)
        self.prt.write('{DASHGO:{N}}'.format(DASHGO=dashgo, N=self.dash_len))

        self.prt.write("{GO_INFO}\n".format(GO_INFO=self.prtfmt.format(**ntgo._asdict())))
        self.gos_printed.add(goid)
        # Do not print hierarchy below this turn if it has already been printed
        if nrp:
            return
        depth += 1
        if self.max_indent is not None and depth > self.max_indent:
            return
        for child in ntobj.children:
            self.prt_hier_rec(child.id, depth)

    @staticmethod
    def _str_dash(depth, single_or_double):
        """Return a string containing dashes (optional) and GO ID."""
        # '-' is default character indicating hierarchy level
        # '=' is used to indicate a hierarchical path printed in detail previously.
        letter = '-' if single_or_double else '='
        return ''.join([letter]*depth)

    def _str_dashgoid(self, ntgo, depth, single_or_double):
        """Return a string containing dashes (optional) and GO ID."""
        dashes = self._str_dash(depth, single_or_double) if self.indent else ""
        return "{DASHES} {GO}{alt:1}".format(DASHES=dashes, GO=ntgo.GO, alt=ntgo.alt)

    def _init_prtfmt(self):
        """Initialize print format."""
        prtfmt = self.gosubdag.prt_attr['fmt']
        prtfmt = prtfmt.replace('{GO} # ', '')
        prtfmt = prtfmt.replace('{D1:5} ', '')
        return prtfmt

#### Examples:
####
#### Print the GO IDs associated with human genes
#### >>> python scripts/wr_hier.py BP --gene2go=gene2go --taxid=9606 --concise -o BP_9606.txt
####
#### Print the hierarchy below Term, GO:0030663
#### >>> python {SCR} GO:0030663
####
#### - GO:0030663	level-05	depth-07	COPI-coated vesicle membrane [cellular_component]
#### -- GO:0012508	level-05	depth-08	Golgi to ER transport vesicle membrane [cellular_component]
#### -- GO:0012509	level-05	depth-08	inter-Golgi transport vesicle membrane [cellular_component]
####
####
#### Write the hierarchy below Term, GO:0030663 into a file
#### >>> python {SCR} GO:0030663 --o=hier_GO_0030663.rpt
####
####   WROTE: hier_GO_0030663.rpt
####
#### Print the hierarchy for biological process, molecular_function, and cellular_component:
#### >>> python {SCR} --o=hier_BP_MF_CC.rpt
####
#### Print hierarchy for BP, MF, CC only printing the first 2 levels.
#### >>> python {SCR} --max_indent=2
#### >>> python {SCR} --max_indent=2 --dash_len=2
####
####
#### Print a conciseened version of the hierarchy for BP, MF, and CC.
#### This will only print a path to a leaf GO Term once.
#### If the path appears a second time, the term is printed again, but its path is not.
#### The presence of a compressed (unprinted) paths is marked by using '=" instead of '-'.
####
####     $ wc -l hier_BP_MF_CC*.rpt
####
####           789583 hier_BP_MF_CC.rpt
####            70152 hier_BP_MF_CC_concise.rpt
####
#### >>> python {SCR} --o=hier_BP_MF_CC_concise.rpt --concise
####
#### Print hierarchy
#### -  26894 GO:0008150	level-00	depth-00	biological_process [biological_process]
#### --    30 GO:0001906	level-01	depth-01	cell killing [biological_process]
#### --   555 GO:0002376	level-01	depth-01	immune system process [biological_process]
#### -- 11208 GO:0065007	level-01	depth-01	biological regulation [biological_process]
####
#### >>> python {SCR}
####
#### This program prints the hierarchy for all GO terms, if no argument is provided.
#### If a GO term is provided as an argument, then the hierarchy of all children
#### for that term is printed.
####
#### """.format(SCR='write_hierarchy')


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.





def write_hier_all(gosubdag, out, root_term):
    """write_hier.py: Prints the entire mini GO hierarchy, with counts of children."""
    out.write('\nTEST ALL: Print all hierarchies:\n')
    objwr = WrHierGO(gosubdag)
    gos_printed = objwr.prt_hier_down(root_term, out)
    print len(gos_printed)
    # assert gos_printed == set(objwr.gosubdag.go2nt)


#################################################################
# Sub-routines to tests
#################################################################
def extract_hier_all(gosubdag, out, root_term, go2geneids):
    """write_hier.py: Prints the entire mini GO hierarchy, with counts of children."""
    # out.write('\nTEST EXTRACTION: Print all hierarchies:\n')
    objwr = WrHierGO(gosubdag,  go2geneids=go2geneids)
    obj = objwr.ext_hier_down(root_term, out)
    return (obj.vertices, obj.edges)



def write_hier_norep(gosubdag, out):
    """Shortens hierarchy report by only printing branches once.
         Prints the 'entire hierarchy' of GO:0000005 the 1st time seen:
           ---     1 GO:0000005    L-02    D-02
           ----     0 GO:0000010   L-03    D-04
         Prints just GO:0000005 (ommit child GO:10) the 2nd time seen:
           ===     1 GO:0000005    L-02    D-02
         '=' is used in hierarchy mark to indicate that the pathes
             below the marked term have already been printed.
    """
    # out.write('\nTEST ALL: Print branches just once:\n')
    objwr = WrHierGO(gosubdag, concise=True)
    gos_printed = objwr.prt_hier_down("GO:0000001", out)
    assert gos_printed == set(objwr.gosubdag.go2nt)


def write_hier_lim(gosubdag, out):
    """Limits hierarchy list to GO Terms specified by user."""
    go_omit = ['GO:0000005', 'GO:0000010']
    go_ids = [go_id for go_id in gosubdag.go2obj if go_id not in go_omit]
    # out.write('\nTEST OMIT: 05 and 10:\n')
    objwr = WrHierGO(gosubdag, include_only=go_ids)
    gos_printed = objwr.prt_hier_down("GO:0000001", out)
    assert not gos_printed.intersection(go_omit), "SHOULD NOT PRINT {GOs}".format(GOs=go_omit)


def write_hier_mrk(gosubdag, out):
    """Print all paths, but mark GO Terms of interest. """
    mark_lst = ['GO:0000001', 'GO:0000003', 'GO:0000006', 'GO:0000008', 'GO:0000009']
    # out.write('\nTEST MARK: 01->03->06->08->09:\n')
    objwr = WrHierGO(gosubdag, go_marks=mark_lst)
    objwr.prt_hier_down("GO:0000001", out)
      #go_marks=[oGO.id for oGO in oGOs_in_cluster])


def fetch_go_hierarcy():
    obo_file_location = os.path.join(constants.GO_DIR, constants.GO_FILE_NAME)
    if not os.path.exists(os.path.join(constants.GO_DIR, constants.GO_FILE_NAME)):
        wget.download(constants.GO_OBO_URL, os.path.join(constants.GO_DIR, constants.GO_FILE_NAME))

    print "Downloading gene-GO associations"
    association_file_location = os.path.join(constants.GO_DIR, constants.GO_ASSOCIATION_FILE_NAME)
    if not os.path.exists(association_file_location):
        wget.download(constants.GO_ASSOCIATION_GENE2GEO_URL,
                      os.path.join(constants.GO_DIR, constants.GO_ASSOCIATION_FILE_NAME))

    print "Loading gene-GO associations"
    # gene2go = download_ncbi_associations(obo_file_location) - why does this line needed?
    go2geneids = read_ncbi_gene2go(association_file_location, taxids=[9606], go2geneids=True)
    geneids2go = read_ncbi_gene2go(association_file_location, taxids=[9606])

    return (go2geneids, geneids2go)

#################################################################
# Driver
#################################################################
def build_hierarcy(roots=['GO:0008150']): #  0008150 0005575 0003674

    go2geneids, geneids2go = fetch_go_hierarcy()

    """Run numerous tests for various reports."""
    dag_fin = os.path.join(constants.GO_DIR, constants.GO_FILE_NAME)
    tic = timeit.default_timer()
    godag = GODag(dag_fin, optional_attrs=['relationship'])
    gosubdag = GoSubDag(godag.keys(), godag)
    toc = timeit.default_timer()
    out = file(os.path.join(constants.BASE_PROFILE, "output", "go_hierarcy.txt"), "w+") # sys.stdout
    dict_result = {}
    for cur_term in roots:
        vertices, edges = extract_hier_all(gosubdag, out, cur_term ,go2geneids)
        # write_hier_norep(gosubdag, out)
        # write_hier_lim(gosubdag, out)
        # write_hier_mrk(gosubdag, out)
        msg = "Elapsed HMS: {}\n\n".format(str(datetime.timedelta(seconds=(toc-tic))))
        sys.stdout.write(msg)
        dict_result[cur_term] = {"vertices" : vertices, "edges": edges}
    return dict_result, go2geneids, geneids2go, get_entrez2ensembl_dictionary()



########################################### ######################
# main
#################################################################
if __name__ == '__main__':
    print build_hierarcy()
