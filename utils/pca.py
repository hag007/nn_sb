
from utils.ensembl2gene_symbol import e2g_convertor
import time
import scipy.special
import matplotlib.pyplot as plt
from matplotlib import style
import matplotlib.ticker as ticker
style.use("ggplot")
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import os
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import proj3d
from matplotlib.lines import Line2D
import random
import matplotlib.cm as cm
import matplotlib.colors as ml_colors
from utils.omic_svm import apply_svm, svm_linear_default, DISTANCE ,make_meshgrid, plot_contours
def plot_pca(dataset, labels_assignment, meta_groups, tested_gene_list_file_name =None):
    actual_labels = list(range(1,len(meta_groups[0])+1))

    labels = [(cur["_name"], cur["_label"]) for i, cur in enumerate(meta_groups[0])]
    labels = [("unknown", 0)]+labels
    X=[]
    y=[]
    for cur in actual_labels:
        X.append(np.average(dataset[np.where(cur == labels_assignment[0]), :][0], axis=0))
        y.append(cur)
    X = np.array(X)
    y = np.array(y)
    pca = PCA(n_components=2)
    labels=[labels[0]] + list(np.array(labels[1:])[~np.isnan(np.average(X, axis=1))])
    y=y[~np.isnan(np.average(X, axis=1))]
    X=X[~np.isnan(np.average(X, axis=1))]
    pca.fit_transform(X[:len(X)].astype(np.float64))

    fig = plt.figure(1, figsize=(20, 20))
    plt.clf()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap='jet')
    for i, x in enumerate(X):
        name = labels[i+1]
        x2, y2, _ = proj3d.proj_transform(x[0], x[1], x[2], ax.get_proj())
        ax.annotate(name,
                    xy=(x2, y2), xytext=(-20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))

    plt.savefig(os.path.join(constants.OUTPUT_GLOBAL_DIR, "PCA_by_average_{}_{}_{}.png").format(constants.CANCER_TYPE,tested_gene_list_file_name.split(".")[0] if tested_gene_list_file_name is not None else "none",time.time()))


def plot_pca_by_samples(dataset, labels_assignment, meta_groups, tested_gene_list_file_name=None, n_components=3):
    actual_labels = list(range(1,len(meta_groups[0])+1))

    X=[]
    y=[]
    for cur in actual_labels:
        for cur_sample in dataset[np.where(cur == labels_assignment[0]), :][0]:
            X.append(cur_sample)
            y.append(cur)
    X = np.array(X)
    y = np.array(y)
    pca = PCA(n_components=n_components)
    pca.fit_transform(X[:len(X)-1].astype(np.float64))

    fig = plt.figure(1, figsize=(20, 20))
    plt.clf()
    if n_components==3:
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap='jet')
    if n_components==2:
        ax = fig.add_subplot(111)
        ax.scatter(X[:, 0], X[:, 1], c=y, cmap='jet')
    plt.savefig(os.path.join(constants.BASE_PROFILE, "output", "PCA_by_samples_{}_{}_{}_{}.png").format(constants.CANCER_TYPE,tested_gene_list_file_name.split(".")[0] if tested_gene_list_file_name is not None else "none", n_components, time.time()))


def plot_detailed_pca(tested_gene_list_file_name, total_gene_list_file_name, gene_expression_file_name,
                      phenotype_file_name, survival_file_name, var_th_index, meta_groups=None,
                      filter_expression=None, feature_names=None, algo_name="", plot_svm=True):

    data = load_integrated_ge_data(tested_gene_list_file_name=tested_gene_list_file_name,
                                   total_gene_list_file_name=total_gene_list_file_name,
                                   gene_expression_file_name=gene_expression_file_name,
                                   phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name,
                                   var_th_index=var_th_index, meta_groups=meta_groups,
                                   filter_expression=filter_expression)
    if data is None:
        print "insufficient data"
        return
    values, h_rows, h_cols, labels_assignment, survival_dataset = data

    X = values
    feature_ids = np.array([x.split('.')[0] for x in h_cols])
    labels = labels_assignment[0]


    n_components = 2
    colormap = cm.jet
    is_transformed=False
    if X.shape[1] < 2:
        print "cannot plot single dimension"
        return

    if X.shape[1] > n_components:
        pca = PCA(n_components=n_components)
        X = pca.fit_transform(X)
        is_transformed=True

    fig = plt.figure(1, figsize=(15, 15))
    plt.clf()
    labels_unique = np.unique(labels)
    myfunc_vec = np.vectorize(lambda x: np.where(labels_unique==x)[0][0])
    colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
                 myfunc_vec(labels) / float(labels_unique.shape[0]-1)]
    if n_components == 3:
        ax = fig.add_subplot(111, projection='3d')
        for x, c, a in zip([x for x in X], colorlist, labels):
            ax.scatter(x[0], x[1], x[2], c=c, cmap='jet', label=a)
    elif n_components == 2:
        ax = fig.add_subplot(111)
        for x, c in zip([x for x in X], colorlist):
            ax.scatter(x[0], x[1], c=c, cmap='jet')
    # colorlist = [ml_colors.rgb2hex(colormap(i)) for i in
    #              np.array(list(range(len(labels)))) / float(len(labels) - 1)]
    colorlist_unique = [ml_colors.rgb2hex(colormap(i)) for i in
                 np.array(list(range(len(labels)))) / float(len(labels_unique) - 1)]
    patches = [Line2D([0], [0], marker='o', color='gray', label=a,
                      markerfacecolor=c) for i, a, c in zip(list(range(labels_unique.shape[0])), labels_unique, colorlist_unique)]
    ax.legend(handles=patches)
    ax.grid(True)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    if is_transformed:
        coeff = np.transpose(pca.components_[0:2, :]) * 150

        top_pca_x_go_terms = feature_ids[np.abs(coeff[:, 0]).argsort()[::-1][:5]]
        top_pca_x_details = [x if feature_names is None else "{}: {}".format(x,feature_names.get(x,"")) for x in top_pca_x_go_terms]
        top_pca_y_go_terms = feature_ids[np.abs(coeff[:, 1]).argsort()[::-1][:5]]
        top_pca_y_details = [x if feature_names is None else "{}: {}".format(x,feature_names.get(x,"")) for x in top_pca_y_go_terms]
        top_pca_go_terms = np.unique(np.append(top_pca_x_go_terms, top_pca_y_go_terms))
        props = dict(boxstyle='round', facecolor='gray', alpha=0.5)
        fig.text(0.01, 0.01,
                 "\n".join(["PC1 features:"] + top_pca_x_details + ["PC2 features:"] + top_pca_y_details),
                 bbox=props)
        n = coeff.shape[0]
        for i in range(n):
            if feature_ids[i] in top_pca_go_terms:
                ax.arrow(0, 0, coeff[i, 0], coeff[i, 1], color='black', alpha=0.5)
                # ax.annotate(feature_ids[i], tuple(coeff[i]), (-100 + random.random() * 200, -100 + random.random() * 200),
                #             textcoords='offset points',
                #             bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                #             arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))
        # for module_i, txt in enumerate(all_hg_score_modules):
        #     ax.annotate(str(txt), (X[:, 0][module_i], X[:, 1][module_i]))

    y=labels

    y = [a - 1 for a in y]
    clf_method = svm_linear_default({'C': [10], 'kernel': ['linear']})
    pr, roc, clf = apply_svm(clf_method, X, y, X, y, DISTANCE)
    if plot_svm:
        xx, yy = make_meshgrid(X[:,0], X[:,1])
        plot_contours(ax, clf, xx, yy,
                      cmap=plt.cm.coolwarm, alpha=0.3)

        ax.set_title('PR score: {}'.format(pr))
        plt.savefig(
            os.path.join(constants.OUTPUT_GLOBAL_DIR, 'pca', constants.DATASET_NAME, algo_name ,"PCA_detailed_{}_{}.png".format(os.path.splitext(tested_gene_list_file_name)[0], time.time())),)

    return (X,labels, pr, roc)
