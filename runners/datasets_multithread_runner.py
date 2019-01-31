import os
from multiprocessing import Process
import utils.aggregate_reports as aggregate_reports
import constants
from runners import bionet_runner
from runners import hotnet2_runner
from runners import jactivemodules_greedy_runner
from runners import jactivemodules_sa_runner
from runners import keypathwayminer_ines_greedy_runner
from runners import netbox_runner
from runners import reactomefi_runner
from runners import matisse_runner
import shutil
from utils.randomize_data import create_random_ds


algo_by_names = {"reactomefi":reactomefi_runner.main,
                 "matisse": matisse_runner.main,
                 "bionet": bionet_runner.main,
                 "keypathwayminer_INES_GREEDY": keypathwayminer_ines_greedy_runner.main,
                 "netbox": netbox_runner.main,
                 "hotnet2": hotnet2_runner.main,
                 "jactivemodules_greedy": jactivemodules_greedy_runner.main,
                 "jactivemodules_sa": jactivemodules_sa_runner.main}

def create_ds_folders(dataset_name):
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "data")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "cache")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "output")))




def run_dataset(dataset_name, expected_genes=None, disease_name=None, score_method=constants.DEG_EDGER, fdr=0.05, algos=None, network_file_name="dip.sif"):
    constants.update_dirs(DATASET_NAME_u=dataset_name)
    if algos is None:
        prcs = [
                Process(target=bionet_runner.main, args=[dataset_name, disease_name, expected_genes, score_method, network_file_name, fdr]),
                Process(target=hotnet2_runner.main, args=[dataset_name, disease_name, expected_genes, score_method, network_file_name]),
                Process(target=netbox_runner.main, args=[dataset_name, disease_name, expected_genes, score_method, network_file_name]),
                Process(target=keypathwayminer_ines_greedy_runner.main, args=[dataset_name, disease_name, expected_genes, score_method, network_file_name]),
                Process(target=reactomefi_runner.main, args=[dataset_name, disease_name, expected_genes, score_method, network_file_name]),
                Process(target=matisse_runner.main, args=[dataset_name, disease_name, expected_genes, score_method, network_file_name])
        ]
    else:
        prcs = []
        for cur in algos:
            if 'jactivemodules' not in cur:
                prcs.append(Process(target=algo_by_names[cur], args=[dataset_name, disease_name, expected_genes, score_method, network_file_name]))
    for cur in prcs:
        cur.start()

    if algos is None or 'jactivemodules_sa' in algos:
        jac_s = Process(target=jactivemodules_sa_runner.main, args=[dataset_name, disease_name, expected_genes, score_method, network_file_name])
        jac_s.start()
        jac_s.join()

    if algos is None or 'jactivemodules_greedy' in algos:
        jac_g = Process(target=jactivemodules_greedy_runner.main, args=[dataset_name, disease_name, expected_genes, score_method, network_file_name])
        jac_g.start()
        jac_g.join()
    for cur in prcs:
        cur.join()





if __name__ == "__main__":
    # datasets=["GE_SOC", "GE_MCF7_2", "GE_TNFa_2", "GE_HC12", "GE_IES", "GE_IEM", "GE_IEN"]
    datasets=["EN_CML", "EN_LICH", "EN_BRCA", "EN_KIRC"] # ["EN_LUNG", "EN_PRAD", "EN_PAAD"]
    algos = ["hotnet2", "bionet", "jactivemodules_greedy", "jactivemodules_sa"]

    for cur_ds in datasets: # datasets: # [1:2]
        print "current folder : {}".format(os.path.basename(cur_ds))
        score_method = constants.PREDEFINED_SCORE
        if cur_ds.startswith("GE"):
            score_method = constants.DEG_EDGER
            if cur_ds.startswith("GE_IE"):
                score_method=constants.DEG_T

        run_dataset(cur_ds, score_method=score_method,
                    algos=algos, network_file_name="dip.sif") #

        # aggregate_reports.aggregate_datasets(os.path.basename(cur_ds))

        # root_path = os.path.join(constants.OUTPUT_GLOBAL_DIR, os.path.basename(cur_ds))
        # print "summary for {}".format(root_path)
        # # #aggregate_reports.aggregate_datasets(os.path.basename(cur_ds))
        # #
        # all_algo_modules = {}
        #
        # for name in os.listdir(root_path):
        #     if os.path.isdir(os.path.join(root_path, name)):
        #      #  modules_report(os.path.basename(cur_ds), name)
        #         modules_summary = pd.read_csv(os.path.join(root_path, name, "modules_summary.tsv"), sep="\t").set_index(
        #             "module")
        #         all_algo_modules[name] = np.array(
        #             modules_summary.index)
        # #
        # emb_modules(os.path.join(constants.OUTPUT_GLOBAL_DIR, os.path.basename(cur_ds)), "all_separated_modules_hg_samples",
        #             "separated_modules_hg_samples", all_algo_modules,
        #             "modules_summary")

#
