import wget
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
        print "Downloading {} to {}".format(link,file_name)
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
    list_of_files_format = ["TCGA.{ds}.sampleMap/HiSeqV2.gz","TCGA.{ds}.sampleMap/{ds}_clinicalMatrix.gz","TCGA.{ds}.sampleMap/miRNA_HiSeq_gene.gz"]
    constants.ALL_CANCER_TYPES = ["PAAD", "OV"] # ["ESCA", "KIRC", "KIRP", "KICH", "LUSC", "LUAD", "COAD", "BRCA", "STAD", "LIHC", "READ", "PRAD", "BLCA", "HNSC", "THCA", "UCEC"]
    for cur in constants.ALL_CANCER_TYPES:
        if cur == "PANCAN": continue
        constants.update_dirs(DATASET_NAME_u="TCGA_"+cur)
        if not os.path.exists(constants.DATA_DIR):
            os.makedirs(constants.DATA_DIR)

        print "fetching data for {} ({}\{})".format(cur,constants.ALL_CANCER_TYPES.index(cur), len(constants.ALL_CANCER_TYPES))
        list_of_files = [fr.format(ds=cur) for fr in list_of_files_format]
        for cur_file_name in list_of_files:
            if not os.path.exists(os.path.join(constants.DATA_DIR,cur_file_name.format(cur))) and not os.path.exists(os.path.join(constants.DATA_DIR,".".join(cur_file_name.split(".")[:-1]))):
                # run_and_printchar(["wget", "https://gdc.xenahubs.net/download/TCGA-{}/Xena_Matrices/{}".format(cur, cur_file_name), constants.TCGA_DATA_DIR])
                download("https://tcga.xenahubs.net/download/"+ cur_file_name, constants.DATA_DIR)
        print "extract data for {}".format(cur)
        for cur_file_name in list_of_files:
            if os.path.exists(os.path.join(constants.DATA_DIR,cur_file_name.split('/')[-1])) and not os.path.exists(os.path.join(constants.DATA_DIR,".".join(cur_file_name.split('/')[-1].split(".")[:-1]))):
                with gzip.open(os.path.join(constants.DATA_DIR,cur_file_name.split('/')[-1]), 'rb') as f_in:
                    with open(os.path.join(constants.DATA_DIR,".".join(cur_file_name.split('/')[-1].split(".")[:-1])), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

        print "delete redundant gz files {}".format(cur)
        for cur_file_name in list_of_files:
            if os.path.exists(os.path.join(constants.DATA_DIR,cur_file_name.split('/')[-1])) and os.path.exists(os.path.join(constants.DATA_DIR,".".join(cur_file_name.split('/')[-1].split(".")[:-1]))):
                os.remove(os.path.join(constants.DATA_DIR,cur_file_name.split('/')[-1]))

        if not os.path.exists(constants.OUTPUT_DIR):
            os.makedirs(constants.OUTPUT_DIR)

        if not os.path.exists(constants.CACHE_DIR):
            os.makedirs(constants.CACHE_DIR)


main()
