from matplotlib import style
style.use("ggplot")
from lifelines.statistics import logrank_test
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import os
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
import time

def km_curve(labels_ids, survival_dataset, tested_gene_expression_headers_columns, gene_group , k=None, label_index=None):
    ax = plt.subplot(111)

    kmf = KaplanMeierFitter()
    all_labels = np.array([y for x in labels_ids for y in x])
    label_event_list = []
    label_duration_list = []
    results = []
    for i, cur_labels in enumerate(labels_ids):
        label_event = survival_dataset[np.in1d(survival_dataset[:, 0], cur_labels) & np.in1d(survival_dataset[:, 0], tested_gene_expression_headers_columns), 4].astype(np.int32)
        label_duration = survival_dataset[np.in1d(survival_dataset[:, 0], cur_labels) & np.in1d(survival_dataset[:, 0], tested_gene_expression_headers_columns), 3].astype(np.int32)
        label_event_list.append(label_event)
        label_duration_list.append(label_duration)
        labels_c = all_labels[~np.in1d(all_labels,cur_labels) & np.in1d(all_labels, tested_gene_expression_headers_columns)]
        label_event_c = survival_dataset[np.in1d(survival_dataset[:, 0], labels_c), 4].astype(np.int32)
        label_duration_c = survival_dataset[np.in1d(survival_dataset[:, 0], labels_c), 3].astype(np.int32)

        lr_results = logrank_test(label_duration, label_duration_c, label_event, label_event_c, alpha=.95)
        if len(label_duration) != 0:
            kmf.fit(list(label_duration), event_observed=list(label_event), label="cluster {} n={}, logrank pval = {}".format(i,len(label_duration), '{0:1.3e}'.format(lr_results.p_value))) # '%.7f' %
            kmf.plot(ax=ax, show_censors=True)
            print "lrank cluster {} vs all: {}".format(i, lr_results.p_value)
            results.append(lr_results.p_value)
            for j, cur_duration in enumerate(label_duration_list[:-1]):
                lr_results = logrank_test(label_duration, label_duration_list[j], label_event, label_event_list[j], alpha=.95)
                print "lrank cluster {} vs cluster {}: {}".format(i, j, lr_results.p_value)
    plt.ylim(0, 1);

    plt.title("clustering survival analysis");
    plt.savefig(os.path.join(constants.BASE_PROFILE,"output" ,"cluster_by_p_{}_{}_k={}_label_i={}_{}.png".format(constants.CANCER_TYPE, gene_group.split("/")[-1],k,label_index , time.time())))
    plt.cla()

    return results

