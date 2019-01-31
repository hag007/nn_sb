import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.io.common import EmptyDataError
from scipy.stats import mannwhitneyu
from scipy.stats import zscore
import constants
import ebm
import matplotlib.colors as ml_colors
from matplotlib.lines import Line2D

hg_report_header = ["GO id", "value", "pval"]
emb_headers= ["algo", "empirical_brown_score", "enriched_terms"]
import seaborn as sns
sns.set(color_codes=True)


def calc_mann_whitney_U(root_path, all_separated_modules_hg_samples_aggregated):
    mannwhitneyu_rank = []
    for name in os.listdir(root_path):
        if os.path.isdir(os.path.join(root_path, name)):
            v_temp = all_separated_modules_hg_samples_aggregated.copy()
            v_temp["empirical_brown_score"] = v_temp["empirical_brown_score"].astype(np.float).apply(
                lambda x: -np.log10(x))
            x = np.array(v_temp[v_temp["algo"] == name]["empirical_brown_score"])
            y = np.array(v_temp[v_temp["algo"] != name]["empirical_brown_score"])
            mannwhitneyu_rank.append({"algo": name, "MWU_pval": mannwhitneyu(x, y).pvalue, "average_score" : x.mean()})
            print "{}".format(mannwhitneyu_rank[-1])
    return mannwhitneyu_rank

def emb_algos(root_path, hg_report_name, general_report_name):

    reports = [os.path.join(root_path, name,hg_report_name)+".tsv" for name in os.listdir(root_path) if os.path.isdir(os.path.join(root_path, name))]
    df_rv = pd.DataFrame() # pd.read_csv(reports[0], sep="\t")[["value"]].rename(columns={'value': reports[0].split("/")[-2]})
    df_pval = pd.DataFrame() # pd.read_csv(reports[0], sep="\t")[["pval"]].rename(columns={'pval': reports[0].split("/")[-2]})
    for hg_report in reports:
        if os.path.isfile(hg_report) and os.path.getsize(hg_report) > 1:
            df_rv = pd.concat([df_rv, pd.read_csv(hg_report, sep="\t")["value"]], axis=1, join='outer').rename(columns={'value': hg_report.split("/")[-2]})
            df_rv = df_rv.fillna(0)
            df_pval = pd.concat([df_pval, pd.read_csv(hg_report, sep="\t")["pval"]], axis=1, join='outer').rename(columns={'pval': hg_report.split("/")[-2]})
            df_pval = df_pval.fillna(1)

    emb_report=[]
    for cur_col in df_pval.columns:
        print "current algo: {}".format(cur_col)
        general_report = pd.read_csv(os.path.join(root_path, cur_col, general_report_name) + ".tsv",sep="\t")
        emb_report.append({
            "algo" : str(cur_col),
            "empirical_brown_score" : str(ebm.KostsMethod(np.array(df_rv, dtype=np.float), np.array(df_pval[cur_col], dtype=np.float))),
            "total_num_genes" : general_report["total_num_genes"][0],
            "enriched_terms" : str((np.sum(df_pval[cur_col]!=1)))})


    pd.DataFrame(emb_report).set_index("algo").to_csv(os.path.join(root_path, hg_report_name+"_aggregated.tsv"), sep="\t")

def emb_modules(root_path, hg_report_name, modules_hg_report_name, algo_to_modules, modules_summary_file_name):

    reports = [os.path.join(root_path, algo,"module_{}_{}.tsv".format(module, modules_hg_report_name)) for algo, modules in algo_to_modules.iteritems() for module in modules]
    df_rv = pd.DataFrame()
    df_pval = pd.DataFrame()
    top_go_terms = np.array([])
    for hg_report in reports:
        new_col_name = hg_report.split("/")[-2] +"_"+ os.path.basename(hg_report).split("_")[1]
        try:
            df_hg_report = pd.read_csv(hg_report, sep="\t")
            top_go_terms=np.append(top_go_terms, df_hg_report[df_hg_report["qval"]<=0.001].index.values)
            # df_hg_report=df_hg_report.loc[df_hg_report["qval"]<=0.001,:]
            df_rv = pd.concat([df_rv, df_hg_report["value"]], axis=1, join='outer').rename(
                columns={'value': new_col_name})
            df_rv = df_rv.fillna(0)
            df_pval = pd.concat([df_pval, df_hg_report["pval"]], axis=1, join='outer').rename(
                columns={'pval': new_col_name})
            df_pval = df_pval.fillna(1)
        except EmptyDataError:
            df_rv[new_col_name] = 0
            df_pval[new_col_name] = 1

    top_go_terms=np.unique(top_go_terms)
    df_rv=df_rv.loc[df_rv.index.isin(top_go_terms),:]
    df_pval=df_pval.loc[df_pval.index.isin(top_go_terms),:]



    emb_report=[]
    print "total number of sig GO terms: {}".format(len(df_pval.index))
    for algo, modules in algo_to_modules.iteritems():
        modules_summary = pd.read_csv(os.path.join(root_path, algo, modules_summary_file_name) + ".tsv", sep="\t").set_index(
            "module")
        for module in modules:
            print "emb algo: {}, module: {}".format(algo, module)
            cur_algo_col = algo+"_"+str(module)
            emb_report.append({
                "algo_module": str(algo)+"_"+str(module),
                "algo" : str(algo),
                "module": str(module),
                "empirical_brown_score" : str(ebm.KostsMethod(np.array(df_rv, dtype=np.float), np.array(df_pval[cur_algo_col], dtype=np.float))),
                "total_num_genes" : modules_summary.loc[[module],:]["#_genes"].iloc[0],
                "enriched_terms" : str((np.sum(df_pval[cur_algo_col]!=1)))})

    pd.DataFrame(emb_report).set_index("algo_module").to_csv(os.path.join(root_path, hg_report_name+"_aggregated.tsv"), sep="\t")



if __name__ == "__main__":




    df_results = pd.DataFrame()
    # datasets = ["GWAS_fasting_insulin", "GWAS_2hr_glucose", "GWAS_adhd", "GWAS_alzheimers", "GWAS_anorexia",
    #             "GWAS_autism", "GWAS_beta-cell_function", "GWAS_bipolar_disorder", "GWAS_blood_pressure_systolic",
    #             "GWAS_body_mass_index", "GWAS_coronary_artery_disease", "GWAS_crohns_disease", "GWAS_cross_disorder"]
    datasets = [name for name in os.listdir(constants.OUTPUT_GLOBAL_DIR) if
                  os.path.isdir(os.path.join(constants.OUTPUT_GLOBAL_DIR, name)) and name.startswith("GWAS_random")] #  and not name.startswith("GWAS_random") and not name.startswith("GWAS_cancer")
    # datasets = ["TNFa_2", "MCF7_2", "SOC", "HC12", "IEM", "IES"]
    # datasets=['GWAS_schizophrenia']
    all_embs = np.array([])

    fig, ax = plt.subplots(figsize=(15, 5))
    for cur_ds in datasets:
        print "current ds: {}".format(cur_ds)
        constants.update_dirs(DATASET_NAME_u=cur_ds)
        root_path = os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME)

        all_algo_modules = {}
        for name in os.listdir(root_path):
            if os.path.isdir(os.path.join(root_path, name)) and name not in  ["data","cache","output"]:
                modules_summary = pd.read_csv(os.path.join(root_path, name, "modules_summary.tsv"), sep="\t")
                if len(modules_summary.index) == 0:
                    continue
                modules_summary=modules_summary.set_index("module")
                all_algo_modules[name] = np.array(
                    modules_summary.index)

        # emb_modules(os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME), "all_separated_modules_hg_samples",
        #             "separated_modules_hg_samples", all_algo_modules,
        #             "modules_summary")
        df_emb = pd.read_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR, constants.DATASET_NAME, "all_separated_modules_hg_samples_aggregated.tsv"), sep='\t')
        all_embs = np.append(all_embs, -np.log10(df_emb['empirical_brown_score'].values))
        # result=calc_mann_whitney_U(root_path, df_emb)


        df_emb = df_emb.iloc[df_emb['empirical_brown_score'].values.argsort(),:]
        df_emb=df_emb.loc[~df_emb['algo'].isin(['matisse', 'reactomefi']),:]
        #algos = np.unique(df_emb['algo'].values)
        algos = np.array(["netbox", "hotnet2", "keypathwayminer_INES_GREEDY", "jactivemodules_sa", "jactivemodules_greedy", "bionet"])
        algos_rank_as_index = df_emb['algo'].apply(lambda x: np.where(algos==x)[0][0]).values
        cmap = plt.cm.jet
        norm = plt.Normalize(vmin=algos_rank_as_index.min(), vmax=algos_rank_as_index.max())
        #

        x = -np.log10(df_emb['empirical_brown_score'].values)
        y = algos_rank_as_index
        s = df_emb['total_num_genes'].values
        ax.scatter(x,y, s, c=[ml_colors.rgb2hex(cmap(a/float(np.size(algos)-1))) for a in y])
        # for a,b,c in zip(x,y,s):
        #     ax.annotate(str(c), (a,b))
        colorlist = [ml_colors.rgb2hex(cmap(a/float(np.size(algos)-1))) for a in np.arange(np.size(algos))]
        patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                          markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
        ax.legend(handles=list(reversed(patches)), loc='upper left', framealpha=0.5)
        ax.set_xlabel("-log10(EMB-pval)", fontdict={"size":20})
        ax.set_ylabel("algos", fontdict={"size":20})
        ax.set_yticks([])
        plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR,cur_ds, "EMB_score_plt.png"))

        # plt.cla()

        # # fig, ax = plt.subplots(figsize=(28, 2))
        # # image = algos_rank_as_index
        # # im= plt.imshow(norm(image.reshape(1,image.shape[0])), aspect='auto', cmap=cmap)
        # # colorlist = [ml_colors.rgb2hex(cmap(i)) for i in norm(np.arange(np.size(algos)))]
        # # patches = [Line2D([0], [0], marker='o', color='gray', label=a,
        # #                   markerfacecolor=c) for i, a, c in zip(list(range(len(algos))), algos, colorlist)]
        # # ax.legend(handles=patches)
        # # ax.grid(False)
        # # plt.show()
        # # plt.imsave(os.path.join(constants.OUTPUT_GLOBAL_DIR,cur_ds,'emb_rank_spectrum.png'), image.reshape(1,image.shape[0]), cmap=cmap)

        # ## ranks = np.zeros(len(result), )
        # # ranks[np.argsort(np.array([cur_res['average_score'] for cur_res in result]))] = np.arange(len(result))[::-1]
        # # for cur_res, cur_rank in zip(result,ranks):
        # #     cur_res['dataset']=cur_ds
        # #     cur_res['rank']=cur_rank
        # #
        # #
        # ## df_results=df_results.append(result)

    plt.clf()
    sns.distplot(all_embs[all_embs> 0],kde=False, bins=50)
    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "emb_dist.png"))
    # plt.cla()
    # # df_results.to_csv(os.path.join(constants.OUTPUT_GLOBAL_DIR,"EMB_MANN.tsv"),sep='\t')


