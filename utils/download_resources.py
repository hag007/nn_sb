import wget
import constants
import os
import gzip
import shutil
import requests
import sys
sys.path.insert(0, '../')
import constants

def download(link, dir_name):
    file_name = os.path.join(dir_name, os.path.basename(link))
    with open(file_name, "wb") as f:
        print "Downloading %s" % file_name
        response = requests.get(link, stream=True)
        total_length = response.headers.get('content-length')

        if total_length is None:  # no content length header
            f.write(response.content)
        else:
            dl = 0
            total_length = int(total_length)
            for data in response.iter_content(chunk_size=4096):
                dl += len(data)
                f.write(data)
                done = int(50 * dl / total_length)
                sys.stdout.write("\r[%s%s]" % ('=' * done, ' ' * (50 - done)))
                sys.stdout.flush()
            sys.stdout.write("\r\n")

def main():
    list_of_files_format = ["TCGA-{}.htseq_counts.tsv.gz","TCGA-{}.htseq_fpkm.tsv.gz","TCGA-{}.htseq_fpkm-uq.tsv.gz","TCGA-{}.GDC_phenotype.tsv.gz","TCGA-{}.survival.tsv.gz","TCGA-{}.mutect2_snv.tsv.gz", "TCGA-{}.mirna.tsv.gz"]
    for cur in constants.ALL_CANCER_TYPES:
        if cur == "PANCAN": continue
        constants.update_dirs(CANCER_TYPE_u=cur)
        if not os.path.exists(constants.TCGA_DATA_DIR):
            os.makedirs(constants.TCGA_DATA_DIR)

        print "fetching data for {} ({}\{})".format(cur,constants.ALL_CANCER_TYPES.index(cur), len(constants.ALL_CANCER_TYPES))
        list_of_files = [fr.format(cur) for fr in list_of_files_format]
        for cur_file_name in list_of_files:
            if not os.path.exists(os.path.join(constants.TCGA_DATA_DIR,cur_file_name.format(cur))) and not os.path.exists(os.path.join(constants.TCGA_DATA_DIR,".".join(cur_file_name.format(cur).split(".")[:-1]))):
                # run_and_printchar(["wget", "https://gdc.xenahubs.net/download/TCGA-{}/Xena_Matrices/{}".format(cur, cur_file_name), constants.TCGA_DATA_DIR])
                download("https://gdc.xenahubs.net/download/TCGA-{}/Xena_Matrices/{}".format(cur, cur_file_name), constants.TCGA_DATA_DIR)
        print "extract data for {}".format(cur)
        for cur_file_name in list_of_files:
            if os.path.exists(os.path.join(constants.TCGA_DATA_DIR,cur_file_name.format(cur))) and not os.path.exists(os.path.join(constants.TCGA_DATA_DIR,".".join(cur_file_name.format(cur).split(".")[:-1]))):
                with gzip.open(os.path.join(constants.TCGA_DATA_DIR,cur_file_name.format(cur)), 'rb') as f_in:
                    with open(os.path.join(constants.TCGA_DATA_DIR,".".join(cur_file_name.format(cur).split(".")[:-1])), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

        print "delete redundant gz files {}".format(cur)
        for cur_file_name in list_of_files:
            if os.path.exists(os.path.join(constants.TCGA_DATA_DIR,cur_file_name.format(cur))) and os.path.exists(os.path.join(constants.TCGA_DATA_DIR,".".join(cur_file_name.format(cur).split(".")[:-1]))):
                os.remove(os.path.join(constants.TCGA_DATA_DIR,cur_file_name.format(cur)))

        if not os.path.exists(constants.OUTPUT_DIR):
            os.makedirs(constants.OUTPUT_DIR)

        if not os.path.exists(constants.CACHE_DIR):
            os.makedirs(constants.CACHE_DIR)


# main()