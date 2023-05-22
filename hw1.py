import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

# function for plotting histograms
def hist(x):
    plt.hist(x, 50)
    plt.show()

# Q1-2-3-4-5
def q2a(count):
    print(len(count))

def q2b(clinical):
    print("(b) How many patient samples are in the dataset? = " + str(len(clinical)))
    print("short: " + str(len(clinical[(clinical == 'short').any(axis=1)])) + ", long: "
          + str(len(clinical[(clinical == 'long').any(axis=1)])))

def q3a(count):
    total = 0
    reads_list = []
    for col in count.columns[1:]:
        total = count[col].sum()
        reads_list.append(total)
    # print(reads_list)
    print(len(reads_list))
    return reads_list

def q3a_graph(reads_list):
    reads = pd.DataFrame(reads_list)
    reads.columns = ['total']
    reads.plot.bar()
    plt.show()

def q3b(count, reads_list):
    n_med = np.median(reads_list)
    tcount_norm = pd.DataFrame(index=count.index, columns=count.columns)
    for col in count.columns[1:]:
        col_idx = count.columns.get_loc(col)-1
        tcount_norm[col] = count[col].apply(lambda x: x * (n_med / reads_list[col_idx]))
    tcount_norm = tcount_norm.iloc[:, 1:]
    # print(tcount_norm)
    # make a new reads list
    reads_list_n = q3a(tcount_norm)
    q3a_graph(reads_list_n)
    return tcount_norm

def q3c(tcount_norm):
    # psuedocount
    tcount_norm = tcount_norm + 1
    count_log = tcount_norm.applymap(np.log2)
    print(count_log)
    # plot histogram
    vals = count_log.stack()
    hist(vals)
    return count_log

def q3d(count_log):
    hist(count_log['TCGA-02-2483-01A-01R-1849-01'])
    hist(count_log['TCGA-02-2485-01A-01R-1849-01'])
    hist(count_log['TCGA-02-2486-01A-01R-1849-01'])
    hist(count_log['TCGA-06-0129-01A-01R-1849-01'])
    hist(count_log['TCGA-06-0178-01A-01R-1849-01'])

def q3e(count_log):
    count_sort = count_log.apply(np.sort, axis=0)
    mean_q = count_sort.mean(axis=1)
    count_qn = pd.DataFrame(index=count_log.index, columns=count_log.columns)
    for sample in count_log.columns:
        sample_data = count_log[sample]
        count_qn[sample] = [mean_q[i] for i in np.searchsorted(count_sort[sample], sample_data)]

    #count_qn = count_sort.apply(lambda x: np.interp(x, x.sort_values(), mean_q))
    print(count_qn)
    # hist(count_qn.stack())
    return count_qn

def q3e_hist(count_qn):
    hist(count_qn['TCGA-02-2483-01A-01R-1849-01'])
    hist(count_qn['TCGA-02-2485-01A-01R-1849-01'])
    hist(count_qn['TCGA-02-2486-01A-01R-1849-01'])
    hist(count_qn['TCGA-06-0129-01A-01R-1849-01'])
    hist(count_qn['TCGA-06-0178-01A-01R-1849-01'])

def q4ab(count_qn):
    short_clinical = clinical[(clinical == 'short').any(axis=1)]
    short_count = pd.DataFrame()
    for col in short_clinical['sampleName']:
        short_count[col] = count_qn[col]
    
    long_clinical = clinical[(clinical == 'long').any(axis=1)]
    long_count = pd.DataFrame()
    for col in long_clinical['sampleName']:
        long_count[col] = count_qn[col]

    p_vals = []
    genes = []
    for index, row in count_qn.iterrows():
        stats, p_val = ranksums(short_count.loc[index], long_count.loc[index])
        p_vals.append(p_val)
    results = pd.DataFrame({'p_value': p_vals})
    sig_genes = results[results['p_value'] < 0.05]
    print(sig_genes.head(10))
    return sig_genes.index

def q4c(sig_genes, deseq2):
    deseq2_index = deseq2[deseq2['pvalue'] < 0.05].index
    print(deseq2_index)
    common_index = sig_genes.intersection(deseq2_index)
    print(common_index) # 197 overlap

def q5a(deseq2):
    num_genes = len(deseq2)
    deseq2['bonferroni'] = deseq2['pvalue'] * num_genes
    sig_genes = deseq2[deseq2['bonferroni'] < 0.05]
    print(sig_genes)
    print(len(sig_genes.index))

def q5b(deseq2):
    p_vals = deseq2['pvalue']
    corrected = multipletests(p_vals, method='fdr_bh')[1]
    deseq2['bh_pvalue'] = corrected
    sig_genes = deseq2[deseq2['bh_pvalue'] < 0.05]
    print(len(sig_genes.index))
    return sig_genes

def q5c(deseq2):
    sorted_genes = deseq2.sort_values('pvalue')
    threshold = 0.05/len(deseq2)
    thresholds = [threshold*(i+1) for i in range(500)]
    fig, ax = plt.subplots()
    ax.plot(range(1, 501), sorted_genes['pvalue'][:500], color='blue', label='DESeq2-p-value')
    ax.plot(range(1, 501), thresholds, linestyle='--', color='red', label='BH-threshold')
    ax.set_xlabel('indices')
    ax.set_ylabel('pvalues')
    ax.set_title('DESeq2 results')
    ax.legend()
    plt.show()

# MAIN
if __name__ == "__main__":
    count = pd.read_csv("HW1-GSE62944-count.csv")
    clinical = pd.read_csv("HW1-GSE62944-clinical.csv")
    deseq2 = pd.read_csv("HW1-DESeq2.csv")
    q2a(count)
    q2b(clinical)
    reads_list = q3a(count)
    q3a_graph(reads_list)
    tcount_norm = q3b(count, reads_list)
    count_log = q3c(tcount_norm)
    q3d(count_log)
    count_qn = q3e(count_log)
    q3e_hist(count_qn)
    sig_genes = q4ab(count_qn)
    q4c(sig_genes, deseq2)
    q5a(deseq2)
    sig_genes = q5b(deseq2)
    q5c(deseq2)
    




