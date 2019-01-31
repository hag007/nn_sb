from infra import *
import constants

constants.update_dirs(CANCER_TYPE_u="SKCM")

tcga = np.array(load_phenotype_data("SKCM_clinicalMatrix"))

gdc = np.array(load_phenotype_data("TCGA-SKCM.GDC_phenotype.tsv"))

old = np.array(load_phenotype_data("SKCM_clinicalMatrix.txt"))

integ = []
integ.append("\t".join(list(gdc[0]) + ["tcga_{}".format(x) for x in tcga[0][1:]] + ["old_{}".format(x) for x in old[0][1:]]))
for cur_gdc in gdc[1:]:
    row = ""
    row += "\t".join(cur_gdc)
    additional = []
    for cur_tcga in tcga[1:]:
        if cur_tcga[0] in cur_gdc[0]:
            additional = cur_tcga[1:]
            break
    additional = "\t".join(additional + ["" for x in range((len(tcga[0][1:]) - len(additional)))])
    row += "\t" + additional
    additional = []
    for cur_tcga in old[1:]:
        if cur_tcga[15] in cur_gdc[0]:
            additional = cur_tcga[1:]
            break
    additional = "\t".join(additional + ["" for x in range((len(tcga[0][1:]) - len(additional)))])
    row +=  "\t" + additional
    integ.append(row)

with open(os.path.join(constants.OUTPUT_DIR, "pheno.txt"),"w+") as f:
    for cur in integ:
        f.write(cur+"\r\n")

