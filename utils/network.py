import pandas as pd
import numpy as np
import os
import constants
import time
import shutil
import sys
import json
import pandas as pd

from df_helpers import to_full_list
from df_helpers import to_full_np
from utils.scripts import format_script
from utils.ensembl2gene_symbol import  e2g_convertor
import zipfile

from utils.go import check_group_enrichment
import go
import multiprocessing
from daemon_multiprocessing import func_star

SH_MODULE_NAME = "module"
SH_NUM_GENES = "#_genes"
SH_ENRICHED = "enriched_groups"
SH_DETAILS = "more_details"
SH_TABLE_HEADERS = [SH_MODULE_NAME, SH_NUM_GENES, SH_ENRICHED, SH_DETAILS]

MODULE_TH = 10

def zipdir(path_to_zip, zip_file_path):
    ziph = zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED)
    for root, dirs, files in os.walk(path_to_zip):
        for file in files:
            ziph.write(os.path.join(root, file))

def get_network_genes(network_name="dip", h_src="ID_interactor_A", h_dst="ID_interactor_B"):
    network_df = pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_name+".sif"), sep="\t")
    src = np.array(network_df[h_src])
    dst = np.array(network_df[h_dst])
    vertices = list(set(np.append(src, dst)))
    return vertices

def remove_subgraph_self_loops(nodes_to_remove, network_file_name=os.path.join(constants.NETWORKS_DIR,"dip.sif"), h_src="ID_interactor_A", h_dst="ID_interactor_B"):
    if len(nodes_to_remove) == 0:
        return network_file_name
    network_df = pd.read_csv(network_file_name, sep="\t")
    filtered_network = network_df[network_df[h_src]!=network_df[h_dst.isin(nodes_to_remove)]]
    new_file_name = os.path.splitext(network_file_name) + "_no_loops" +".sif"
    filtered_network.to_csv(new_file_name, sep="\t", index=False)
    return filtered_network

def remove_subgraph_by_nodes(nodes_to_remove, network_file_name=os.path.join(constants.NETWORKS_DIR,"dip.sif"), h_src="ID_interactor_A", h_dst="ID_interactor_B", ts=str(time.time())):
    if len(nodes_to_remove) == 0:
        return network_file_name
    network_df = pd.read_csv(network_file_name, sep="\t")
    filtered_network = network_df[~(network_df[h_src].isin(nodes_to_remove) | network_df[h_dst].isin(nodes_to_remove))]
    new_file_name = os.path.splitext(network_file_name)[0] + ts +".sif"
    filtered_network.to_csv(new_file_name, sep="\t", index=False)
    return new_file_name



def summary_intergrative_reports(all_hg_reports, modules_summary, total_hg_report, algo_name, module_genes, disease_name, expected_genes, report_file_name, dataset_name):



    general_algo_report(algo_name, all_hg_reports, module_genes, modules_summary, report_file_name, total_hg_report, dataset_name)

    if constants.DISEASE_MODE:
        disease_algo_report(algo_name, disease_name, expected_genes, module_genes, modules_summary, report_file_name)

    if constants.EMB_MODE:
        emb_score_report(algo_name, report_file_name, "hg_samples", total_hg_report, dataset_name)


def emb_score_report(algo_name, report_file_name, hg_sample_file_name, hg_report, dataset_name):
    samples = [{go.HG_GO_ID : cur_term[go.HG_GO_ID], go.HG_GO_NAME : cur_term[go.HG_GO_NAME], go.HG_GO_ROOT : cur_term[go.HG_GO_ROOT], go.HG_VALUE : cur_term[go.HG_VALUE], go.HG_PVAL : cur_term[go.HG_PVAL], go.HG_QVAL : cur_term[go.HG_QVAL] } for cur_term in hg_report]
    df_emb = pd.DataFrame(samples)
    df_emb.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, dataset_name, algo_name,
                     "{}_{}.tsv".format(report_file_name, hg_sample_file_name)),sep="\t", index=False)
    return df_emb


def disease_algo_report(algo_name, disease_name, expected_genes, module_genes, modules_summary, report_file_name, dataset_name):

    disease_data = {
        "disease_name": disease_name,
        "num_of_modules": len(modules_summary),
        "TP+FN_(_true_)": len(expected_genes),
        "TP+TN_(_retrieved_)": len(module_genes),
        "TP/(TP+TN)_(_precision_)": 0,
        "TP/(TP+FN)_(_recall_)": 0,
        "F1": 0,
        "TP": 0,
        "module_size_avg" : 0,
        "module_size_std" :0
    }
    if len(modules_summary) > 0:
        modules_summary = pd.DataFrame(modules_summary)
        disease_genes_extracted = float(len(set(module_genes).intersection(expected_genes)))
        disease_data["TP"] = disease_genes_extracted
        disease_data["TP/(TP+TN)_(_precision_)"] = disease_genes_extracted / len(module_genes)
        disease_data["TP/(TP+FN)_(_recall_)"] = disease_genes_extracted / len(expected_genes)
        if (disease_data["TP/(TP+TN)_(_precision_)"] + disease_data["TP/(TP+FN)_(_recall_)"]) == 0:
            disease_data["F1"] = 0
        else:
            disease_data["F1"] = 2 * ((disease_data["TP/(TP+TN)_(_precision_)"] * disease_data["TP/(TP+FN)_(_recall_)"]) /
                                  (disease_data["TP/(TP+TN)_(_precision_)"] + disease_data["TP/(TP+FN)_(_recall_)"]))

        disease_data["module_size_avg"] = modules_summary[SH_NUM_GENES].mean()
        disease_data["module_size_std"] = modules_summary[SH_NUM_GENES].std()


    pd.DataFrame([disease_data]).to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, dataset_name, algo_name, "{}_disease.tsv".format(report_file_name)),sep="\t", index=False)


def general_algo_report(algo_name, all_hg_reports, module_genes, modules_summary, report_file_name, total_hg_report, dataset_name):
    data = {}
    if len(modules_summary) > 0 :
        df_summary = pd.DataFrame(modules_summary)
        data = {"num_of_modules": df_summary.index.size,
                "module_size_avg": df_summary[SH_NUM_GENES].mean(),
                "module_size_std": df_summary[SH_NUM_GENES].std(),
                "total_num_genes": len(module_genes)
                }


    if len(all_hg_reports) > 0:
        df_all_hg = [pd.DataFrame(x) for x in all_hg_reports]
        enrichment_dist = [x.index.size for x in df_all_hg]
        pval_dist = [np.array(x[go.HG_PVAL].astype(np.float).apply(lambda x: -np.log10(x))) if x.index.size > 0 else np.array([]) for x in df_all_hg]
        modules_enrichment_data = {"module_enriched_terms_avg": np.average(enrichment_dist),
                                   "module_enriched_terms_std": np.std(enrichment_dist),
                                   "module_enriched_terms_signal_avg_avg": np.average([np.average(x) if len(x) > 1 else 0
                                                                                        for x in pval_dist]),
                                   "module_enriched_terms_signal_avg_std": np.average([np.std(x) if len(x) > 1 else 0
                                        for x in pval_dist]),
                                   "module_enriched_terms_signal_std_avg": np.std([np.average(x) if len(x) > 1 else 0
                                        for x in pval_dist]),
                                   "module_enriched_terms_signal_std_std": np.std([np.std(x) if len(x) > 1 else 0
                                        for x in pval_dist])}


        data.update(modules_enrichment_data)
        data["module_enriched_terms_signal_score"] = \
            data['module_enriched_terms_signal_avg_avg'] / ((data[
                                                                 'module_enriched_terms_signal_avg_std'] +
                                                             data[
                                                                 'module_enriched_terms_signal_std_avg'] +
                                                             data[
                                                                 'module_enriched_terms_signal_std_std']) * \
                                                            (data[
                                                                 "total_num_genes"] /
                                                             data[
                                                                 "module_size_avg"] *
                                                             data[
                                                                 "num_of_modules"]))

    if len(total_hg_report) > 0:
        df_total_hg = pd.DataFrame(total_hg_report)
        all_enrichment_data = {
            "total_enriched_terms_avg" : df_total_hg[go.HG_PVAL].astype(np.float).apply(lambda x: -np.log10(x)).mean(),
            "total_enriched_terms_std" : df_total_hg[go.HG_PVAL].astype(np.float).apply(lambda x: -np.log10(x)).std(),
            "total_num_enriched_terms": len(total_hg_report)
        }
        data.update(all_enrichment_data)

    df = pd.DataFrame()
    if len(data) >0:
        df = pd.DataFrame([data])

    df.to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, dataset_name, algo_name,
                     "{}_general.tsv".format(report_file_name)), sep="\t", index=False)



def output_modules(output_file_name, modules, score_file_name, output_base_dir=""):
    output_data = create_modules_output(modules, score_file_name)
    file(output_file_name, 'w+').write(output_base_dir + "\n")
    json.dump(output_data, file(output_file_name, 'a+'))
    sys.stdout.write(output_file_name)

def reduce_to_dict(x,y):
    if y["id"] in x:
        x[y["id"]]["modules"] = x[y["id"]]["modules"] + y["modules"]
    else:
        x[y["id"]]=y
    return x

def merge_two_dicts(x, y):

    z = x.copy()
    z.update(y)
    return z

def create_modules_output(modules, score_file_name):
    scores=None
    if score_file_name is not None:
        scores = pd.read_csv(score_file_name,sep="\t").set_index("id")

        if constants.IS_PVAL_SCORES:
            scores["score"] = scores["pval"].apply(lambda x: -np.log10(x))

    zero_scores = [ {"score" : 0, "id" : gene} for module in modules for gene in module if scores is None or gene not in scores.index]
    if len(zero_scores) !=0:
        zero_scores = pd.DataFrame(zero_scores).set_index("id")
        zero_scores=zero_scores[~zero_scores.index.duplicated(keep='first')]
        scores = pd.concat([scores, zero_scores],axis=0)
    return [merge_two_dicts({"id" : k}, v) for k,v in reduce(reduce_to_dict, [{"eid": gene, "modules": [i], "id": gene, "gene_symbol": e2g_convertor([gene])[0], "score" : scores.loc[gene,"score"]} for i, module in enumerate(modules) for gene in module],\
            {}).iteritems()]

def draw_network(modules, score_file_name, network_file_name, h_src="ID_interactor_A", h_dst="ID_interactor_B"):
    active_genes = [y for x in modules for y in x]
    output = [{"data" : x, "label" : x["eid"], "selected" : True } for x in create_modules_output(modules, score_file_name)]
    active_edges = [[x[h_src], x[h_dst]] for i, x in pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_file_name), sep="\t").iterrows() if x[h_src] in active_genes and x[h_dst] in active_genes]
    additional_edges = [[x[h_src], x[h_dst]] for i, x in pd.read_csv(os.path.join(constants.NETWORKS_DIR, network_file_name), sep="\t").iterrows() if not (x[h_src] in active_genes and x[h_dst] in active_genes) and (x[h_src] in active_genes or x[h_dst] in active_genes)]
    additional_nodes = [y for x in (active_edges + additional_edges) for y in x if y if y not in active_genes]
    additional_nodes = list(set(additional_nodes))

    return output + [{"data" : {"id" : x, "eid" : x, "modules" : []}, "label" : ""} for x in additional_nodes] + [{"data": {"id" : x[0]+"_"+x[1], "source":x[0], "target":x[1]}, "label" : ""} for x in additional_edges] + [{"data": {"id" : x[0]+"_"+x[1], "source":x[0], "target":x[1]}, "label" : "-"} for x in active_edges]



def generate_report_from_template(output_file_name, cy, algo_name="", hg_report=[], dataset_name="" ,
                                  disease_genes_statistics=[], modules_summary=[]):

    hg_report = to_full_list(hg_report, "#")
    disease_genes_statistics = to_full_list(disease_genes_statistics, "#")
    modules_summary = to_full_list(modules_summary, "#")

    report_file_name=format_script(os.path.join(constants.TEMPLATES_DIR, "graph.html"),
                  DATA=json.dumps(cy), HG_REPORT=json.dumps(hg_report),
                  MODULES_SUMMARY=json.dumps(modules_summary), NUM_OF_GENES=len([x for x in cy if not x["data"].has_key("source") and len(x["data"]["modules"])>0]),
                  DISEASE_GENES=json.dumps(disease_genes_statistics))
    output_dir = os.path.join(constants.OUTPUT_GLOBAL_DIR, dataset_name, algo_name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    shutil.move(report_file_name,
                os.path.join(output_dir, "graph_{}.html".format(output_file_name)))
    return "graph_{}.html".format(output_file_name)






def build_all_reports(algo_name, dataset_name, modules, all_bg_genes, score_file_name, network_file_name, disease_name=None, expected_genes=None):

    output_base_dir = os.path.join(constants.OUTPUT_GLOBAL_DIR, dataset_name, algo_name)
    if not os.path.exists(output_base_dir):
        os.makedirs(output_base_dir)

    manager=multiprocessing.Manager()
    all_hg_reports = manager.list()
    modules_summary = manager.list()

    params=[]
    p=multiprocessing.Pool(5)
    for i, module in enumerate(modules):
        params.append([module_report, [algo_name, i, module, all_bg_genes[i], score_file_name, network_file_name, dataset_name, all_hg_reports,
             modules_summary]])
        # module_report_star([algo_name, i, module, all_bg_genes[i], score_file_name, network_file_name, dataset_name, all_hg_reports, modules_summary])

    p.map(func_star, params)

    modules_summary=list(modules_summary)
    all_hg_reports=list(all_hg_reports)

    # module_genes = list(set(reduce((lambda x, y: x + y), modules, [])))
    modules_larger_than_k, module_larger_than_k_genes, k_hg_reports, k_modules_summary = \
        get_k_threshold_modules(modules, all_hg_reports, modules_summary)

    df_summary = pd.DataFrame([], columns=['#_genes'])
    df_summary.index.name="module"
    bg_genes = []
    if len(modules) > 0:
        df_summary=pd.DataFrame(modules_summary).set_index("module")
        bg_genes = all_bg_genes[0]

    df_summary.to_csv(
        os.path.join(constants.OUTPUT_GLOBAL_DIR, dataset_name, algo_name, "modules_summary.tsv"), sep="\t")
    generate_algo_report(algo_name, modules, bg_genes, all_hg_reports, disease_name, expected_genes,
                         modules_summary, score_file_name, network_file_name, "all_modules" ,dataset_name)

    generate_algo_report(algo_name, modules_larger_than_k, bg_genes, k_hg_reports, disease_name, expected_genes,
                         k_modules_summary, score_file_name, network_file_name, "k_{}_modules".format(MODULE_TH) ,dataset_name)

    return output_base_dir



def get_k_threshold_modules(modules, all_hg_reports, modules_summary):
    modules_larger_than_k = [cur for cur in modules if len(cur) >= MODULE_TH]
    module_larger_than_k_genes = list(set(reduce((lambda x, y: x + y), modules_larger_than_k, [])))
    k_modules_summary = [modules_summary[i] for i, cur in enumerate(modules) if len(cur) >= MODULE_TH]
    k_hg_reports = []
    if constants.HG_MODE:
        k_hg_reports = [all_hg_reports[i] for i, cur in enumerate(modules) if len(cur) >= MODULE_TH]
    return modules_larger_than_k, module_larger_than_k_genes, k_hg_reports, k_modules_summary


def generate_algo_report(algo_name, modules, bg_genes, all_hg_reports, disease_name, expected_genes,
                         modules_summary, score_file_name, network_file_name, report_name, dataset_name):
    hg_report = []
    module_genes = list(set([gene for module in modules for gene in module]))
    if (constants.HG_MODE or constants.EMB_MODE) and constants.ALGO_HG_MODE:
        hg_report = check_group_enrichment(module_genes, bg_genes, algo_name)
    cy = draw_network(modules, score_file_name, network_file_name)
    generate_report_from_template(report_name, cy, algo_name, hg_report, dataset_name, [], modules_summary)
    summary_intergrative_reports(all_hg_reports, modules_summary, hg_report, algo_name, module_genes, disease_name,
                                 expected_genes, report_name, dataset_name)


def module_report(algo_name, module_index, module, bg_genes, score_file_name, network_file_name, dataset_name, all_hg_reports=None, modules_summary=None):
    print "summarize module {} for algo {} and dataset {}".format(module_index, algo_name, dataset_name)

    file(os.path.join(constants.OUTPUT_DIR, "{}_module_genes_{}.txt".format(algo_name, module_index)), "w+").write(
        "\n".join(module))
    file(os.path.join(constants.OUTPUT_DIR, "{}_bg_genes_{}.txt".format(algo_name, module_index)), "w+").write(
        "\n".join(bg_genes))
    modules_summary_row = {SH_MODULE_NAME: module_index, SH_NUM_GENES: len(module)}
    hg_report = []
    if constants.HG_MODE:
        hg_report = check_group_enrichment(list(module), list(bg_genes), algo_name, str(module_index))
        modules_summary_row[SH_ENRICHED] = len(hg_report)
        emb_score_report(algo_name, "module_" + str(module_index), "separated_modules_hg_samples", hg_report, dataset_name)
    cy = draw_network([[] for a in range(module_index)] + [module], score_file_name, network_file_name)
    report_output_file_name = generate_report_from_template(algo_name + str(module_index), cy, algo_name, hg_report, dataset_name)
    modules_summary_row[SH_DETAILS] = report_output_file_name
    if all_hg_reports is not None and modules_summary is not None:
        all_hg_reports.append(hg_report)
        modules_summary.append(modules_summary_row)
    return hg_report, modules_summary_row


