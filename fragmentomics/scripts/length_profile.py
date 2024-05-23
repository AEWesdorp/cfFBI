import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['font.family'] = ['Arial']

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def get_normalized_smoothing(cyc_dict, sample, smoothing_factor =1, read_length_range=(31, 1000)):
    d = defaultdict(int)
    for key in range(read_length_range[0], read_length_range[1]):
        # if smoothing_factor == 1, no smoothing is performed.
        written_key = key // smoothing_factor * smoothing_factor
        # If there is no value with associated length, add 0 to the key.
        if key % smoothing_factor != 0:
            try:
                d[written_key] += cyc_dict['amount'][key]
            except KeyError:
                d[written_key] += 0
        else:
            try:
                d[written_key] += cyc_dict['amount'][key]
            except KeyError:
                d[written_key] += 0
    rl = pd.DataFrame(d, index = ['amount']).T
    rl[sample] = rl['amount'].apply(lambda x: x / rl['amount'].sum())
    return rl


def get_combined_normalized_read_length(path_csv,
                                        suffix,
                                        samples,
                                        smoothing_factor = 1,
                                        read_length_range: tuple=(31, 1000)):
    all_rl = []
    for sample in samples:
        cfdna =pd.read_csv(f"{path_csv}/{sample}{suffix}",sep = "\t",comment = "#", skiprows=11, names = ['len', 'amount', 'inward','outward','other'] )
        cfdna = cfdna[['len','amount']]
        cfdna.index = cfdna['len']
        cfdna_dict = cfdna.to_dict()
        assert smoothing_factor >= 1
        rl = get_normalized_smoothing(cfdna_dict, sample,  smoothing_factor=smoothing_factor, read_length_range=read_length_range)
        all_rl.append(rl)
    df = pd.concat(all_rl, axis = 1)
    plasma_set = df[samples].T
    return plasma_set


def get_colors_list(samples):
    color_list = []
    color_dict = {'119':'black','210':'red','240':'blue','OES':'#EC9D40','EAC':'#EC9D40', '10percent':'black', '1percent':'black', '1':'black',  '2percent':'black','05percent':'black'}
    for sample in samples:
        color_list.append(color_dict[sample[0:3]])
    return color_list


def plot_length(ps, samples, out_path, figsize=(3, 8), colors=False):

    ## Plot components:
    if colors:
        pass
    else:
        colors = get_colors_list(samples)
    i = 0
    fig, axes = plt.subplots(len(samples), figsize=figsize, sharex=True)
    #axes[0].axvline(167, color='black',linewidth=0.5)

    for sample, color in zip(samples, colors):
        axes[i].plot(ps.keys(), ps.loc[sample], label=sample, color=color)
        #     axes[i].legend(loc=(1.01,0),frameon=False)
        axes[i].spines['top'].set_visible(False)
        #     axes[i].spines['bottom'].set_visible(False)
        #axes[i].spines['left'].set_visible(False)
        axes[i].spines['right'].set_visible(False)
        #     axes[i].set_ylim(0, 0.016)
        #axes[i].set_yticks([])
        axes[i].set_ylabel(sample, rotation=0)
        # axes[i].s('')
        i += 1
    axes[i - 1].set_xlabel('cfDNA length')
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, transparent=True)
    plt.savefig(out_path + ".png", dpi=300, transparent=True)
    plt.show()


if __name__ == '__main__':
    out_path = '/Users/lchen/00_projects/FBI/output/04_fragmentomics/MT.pdf'
    path_csv = '/Users/lchen/00_projects/FBI/output/02_tables/00_pipeline_postprocess/read_length_stats/MT_read_length'
    #suffix = "_alligned_host_pe_srt_EquCab3_major_chromosomes.txt_readlength.txt"
    suffix =  '_alligned_host_pe_srt_EquCab3_MT.txt_readlength.txt'
    samples = ['119089','119352','2104671','2105097','2401142']
    df_CYC = get_combined_normalized_read_length(path_csv,
                                             suffix,
                                             samples,
                                             smoothing_factor = 1,
                                             read_length_range=(31,1600) )
    ## plot:
    plot_length(df_CYC.iloc[:, :200], samples, out_path, figsize = (8,15))


