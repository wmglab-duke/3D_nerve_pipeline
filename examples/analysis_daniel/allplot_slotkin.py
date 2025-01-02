import json
import os

import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import stats
from scipy.stats import pearsonr

os.chdir('../../')
import sys

sys.path.pop(-2)
from src.core.plotter import datamatch_merge

mpl.rcParams['figure.dpi'] = 400
sns.set_style('whitegrid')
pd.options.mode.chained_assignment = None


def corrcalc(data, comparison):
    corrs = []
    for nsim in pd.unique(data['nsim']):
        # ax.set_title(f'fiber diam: {s}μm')
        corr = {}
        for sample in pd.unique(data['sample']):
            thisdat = data[(data["nsim"] == nsim) & (data["sample"] == sample)]
            corr[sample] = round(pearsonr(thisdat[comparison[0]], thisdat[comparison[1]])[0], 3)
        corrs.append(corr)
    return corrs


def addpwfd(data, sim):
    with open('examples/analysis/plotconfig.json') as f:
        config = json.load(f)
    nsim_key = config['sim_data'][sim]['nsim_key']
    for nsim in nsim_key:
        pulse_width = nsim_key[nsim]['pulse_width']
        fiber_diam = nsim_key[nsim]['fiber_diam']
        nsim = int(nsim)
        data.loc[data.nsim == nsim, 'pulse_width'] = pulse_width
        data.loc[data.nsim == nsim, 'fiber_diam'] = fiber_diam
    return data


def pe(correct, est):
    """Calculate the percent error.

    :param correct: correct value
    :param est: estimated value
    :return: percent error
    """
    return 100 * abs(est - correct) / correct


def pe_noabs(correct, est):
    """Calculate the percent error.

    :param correct: correct value
    :param est: estimated value
    :return: percent error
    """
    return 100 * (est - correct) / correct


import math


def compute_reorder_cost(order1, order2):
    assert set(order1) == set(order2)
    assert len(set(order1)) == len(order1)
    assert len(set(order2)) == len(order2)
    sumdist = 0
    for item in order1:
        distance = abs(order1.index(item) - order2.index(item))
        sumdist += distance
    if len(order1) % 2 == 0:
        maxlen = sum(range(len(order1))[int(math.ceil(0.5 * len(order1))) :]) * 2 / len(order1)
    else:
        maxlen = (
            range(len(order1))[int(math.floor(0.5 * len(order1)))]
            + (sum(range(len(order1))[int(math.ceil(0.5 * len(order1))) :])) * 2
        ) / len(order1)
    return sumdist / len(order1) / maxlen


def calculate_dose_response(df, threshold_column, outcolumn, grouping_columns):
    def calculate_dose_response_per_sample(group):
        threshold_values = group[threshold_column].values
        activated_fibers = group[threshold_column].values[:, None] <= threshold_values
        dose_response = activated_fibers.sum(axis=0) / len(group)
        group[outcolumn] = dose_response
        return group

    return df.groupby(grouping_columns).apply(calculate_dose_response_per_sample)


def rownames(g, row_template, **ylabel_kws):
    args = dict(row_var=g._row_var)
    # Draw the row titles on the left edge of the grid
    for i, row_name in enumerate(g.row_names):
        ax = g.axes[i, 0]
        args.update(dict(row_name=row_name))
        title = row_template.format(**args)
        ax.set_ylabel(title, **ylabel_kws)


gogo = "initial"

cath_comparison = ["cathodic", '3D', 3]
center_comparison = ["center", '3D', 3]
an_comparison = ["anodic", '3D', 3]
main_comparison = center_comparison

simNUM = input('Gimme simNUM: ')
if simNUM == 10:
    main_comparison = cath_comparison

pal2d3d = ['#d95f02', '#7570b3']

# code whole file to optionally run on 3 or 10 so that I can run all monpolar data
# %% set which comparison to run
print('threshload')
# base data
threshload = pd.read_csv(f'thresh_unmatched_sim{simNUM}_og.csv')

center = threshload['sample'].astype(str).str[2] == '1'
threshload.loc[center, 'contact'] = 'center'
threedcontact = threshload['sample'].astype(str).str[2] == '3'
threshload.loc[threedcontact, 'contact'] = '3D'
# set all inners and outers where 'type' is 3D to 0
threshload.loc[(threshload['type'] == '3D'), 'inner'] = 0
threshload.loc[(threshload['type'] == '3D'), 'outer'] = 0
# %% inners need to match the cathodic leading contact
# Query the rows with type '3D'
df_3d = threshload.query("type == '3D'")

# Query the rows with type '2D' and sample ending in '1'
df_2d = threshload.query(f"type == '2D' and contact in {main_comparison}")

# Merge the 3D and 2D data, keeping track of original row indices
merged_df = pd.merge(
    df_3d, df_2d, on=['nerve_label', 'master_fiber_index', 'nsim'], suffixes=('_3d', '_2d'), how='left'
)  # TODO remove this how

# Update the 'inner', 'outer', and 'fiber' columns in the original DataFrame
threshload.loc[df_3d.index, 'inner'] = merged_df['inner_2d'].values
threshload.loc[df_3d.index, 'outer'] = merged_df['outer_2d'].values
threshload.loc[df_3d.index, 'fiber'] = merged_df['fiber_2d'].values
# %%
if gogo == "initial":  # remove all where nerve_label length is >2
    threshload = threshload[threshload['nerve_label'].str.len() < 3]
threshload['type'] = threshload['type'].replace({'2D': '2DEM', '3D': '3DM'})
threshload = addpwfd(threshload, str(simNUM))
threshload['fiber_diam'] = threshload['fiber_diam'].astype(int)
# elif gogo=="deformedasc": #remove all where nerve_label does not contain "asc"
#     threshload = threshload[threshload['nerve_label'].str.contains("asc")]
#     #check the third digit of the sample number is 2, in that case, "contact" is cathodic. If 0, anodic
threshload['contact'] = (
    threshload['sample'].astype(str).str[2].replace({'0': 'anodic', '2': 'cathodic', '3': '3D', '1': 'center'})
)
threshload.sort_values(by='sample')

# %%Setup
print('matched')
# 3DM and 2DEM trhresholds as different col same row
matched = datamatch_merge(
    threshload.query('type=="2DEM"'),
    threshload.query('type=="3DM"'),
    'threshold',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='type')
# for matched boxplots
repeated = matched.melt(
    id_vars=[x for x in matched.columns if 'threshold' not in x],
    value_vars=['threshold3d', 'threshold'],
    var_name='type',
    value_name='threshold',
)
repeated['type'] = repeated['type'].replace({'threshold': '2DEM', 'threshold3d': '3DM'})
repeated['type'] = pd.Categorical(repeated['type'], categories=['2DEM', '3DM'], ordered=True)
# %%
print('dose-response')
# dose response
drdat = threshload.copy().rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
drdat = calculate_dose_response(
    drdat,
    'threshold',
    'percent_activated',
    grouping_columns=['modeltype', 'samplenum', 'fiber_diam', 'sim'],
)
drdat.sort_values('modeltype', inplace=True)
drmatch = datamatch_merge(
    drdat.query('modeltype=="2DEM"'),
    drdat.query('modeltype=="3DM"'),
    'percent_activated',
    merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
).drop(columns='modeltype')

# %% BEGIN DEFORMATION ANALYSIS
print('deformation')
# match fascicle thresholds
sns.set(font_scale=1.75)
sns.set_style('whitegrid')
threshes = pd.concat(
    [pd.read_csv(f'thresh_unmatched_sim{simNUM}_og.csv'), pd.read_csv(f'thresh_unmatched_sim{simNUM}_def.csv')]
)
# remove all rows where nerve label contains "asc" and "sample" contains "3"
threshes = threshes[~((threshes['nerve_label'].str.contains('asc')) & (threshes['sample'].astype(str).str[2] == '3'))]
# duplicate all samples where nerve label contains "def" and "sample" contains "3"
defdat = threshes[(threshes['nerve_label'].str.contains('def')) & (threshes['sample'].astype(str).str[2] == '3')]
# replacee all "1" with "9" in sample column
defdat['sample'] = defdat['sample'].astype(str).str.replace('1', '9').astype(int)
# replace all "def" with "asc" in nerve label column
defdat['nerve_label'] = defdat['nerve_label'].str.replace('def', 'asc')
# combine with original data
threshes = pd.concat([threshes, defdat])
threshes['deformation'] = None
# where nerve label contains "def" deformation is "3D-3D"
threshes.loc[threshes['nerve_label'].str.contains('def'), 'deformation'] = '3D-3D'
# where nerve label contains "asc" deformation is "2D-3D"
threshes.loc[threshes['nerve_label'].str.contains('asc'), 'deformation'] = '2D-3D'
# else deformation is "none"
threshes.loc[threshes['deformation'].isna(), 'deformation'] = "Undeformed"
# strip "def" and "asc" from nerve labels
threshes['nerve_label'] = threshes['nerve_label'].str.replace('def', '').str.replace('asc', '')

newdefdat = threshes.copy()

# remove all nerve_label = '6L'
newdefdat = newdefdat[newdefdat['nerve_label'] != '6L']
newdefdat = newdefdat[newdefdat['nerve_label'] != '2R']

newdefdat['type'] = newdefdat['type'].replace({'2D': '2DEM', '3D': '3DM'})
newdefdat = addpwfd(newdefdat, simNUM)
newdefdat['fiber_diam'] = newdefdat['fiber_diam'].astype(int)
# contact is cathodic if the third digit of sample int is 2, anodic if 0
newdefdat['contact'] = (
    newdefdat['sample'].astype(str).str[2].replace({'2': 'cathodic', '1': 'center', '0': 'anodic', '3': '3D'})
)

# set deformation as ordered categorical
newdefdat['deformation'] = pd.Categorical(
    newdefdat['deformation'], categories=['Undeformed', '2D-3D', '3D-3D'], ordered=True
)
newdefdat['nerve_label'] = pd.Categorical(newdefdat['nerve_label'], categories=['2L', '3R', '5R', '6R'], ordered=True)
newdefdat['deformed'] = newdefdat['deformation'] != "Undeformed"
deftomatch = newdefdat.copy()

# remove unused colors from palette
defpal = [sns.color_palette('colorblind')[ind] for ind in [0, 2, 3, 5]]
defdefcomp = newdefdat.query('type=="3DM" and deformation != "2D-3D"')
defdefcomp['deformed'] = defdefcomp['deformation'] != "Undeformed"
defsamples = ["2L", "3R", "5R", "6R"]
# %% new idea analysis
concats = []
for deftype in ["Undeformed", "3D-3D"]:
    thismatch = deftomatch.query('deformation==@deftype')
    matchednow = datamatch_merge(
        thismatch.query('type=="2DEM"'),
        thismatch.query('type=="3DM"'),
        'threshold',
        merge_cols=["model", "sim", "nerve_label", "nsim", "master_fiber_index"],
    ).drop(columns='type')
    matchednow['deformation'] = deftype
    concats.append(matchednow)
concats = pd.concat(concats)
concats['deformation'] = pd.Categorical(concats['deformation'], categories=['Undeformed', '3D-3D'], ordered=True)
concats.to_csv(r'C:\Users\dpm42\My Drive\ClassFilesGrad\SLOTKIN\analysis\matched.csv', index=False)

# %% dose-response
print('deformdr')
defdr = newdefdat.copy().rename(columns={'sample': 'samplenum', 'type': 'modeltype'})
defdr = calculate_dose_response(
    defdr,
    'threshold',
    'percent_activated',
    grouping_columns=['modeltype', 'samplenum', 'fiber_diam', 'sim', 'deformation', 'pulse_width'],
)
sys.exit('prepdone')
# %% activation order
newdefdr = defdr.copy()
newdefdr = newdefdr.query("deformation=='3D-3D'").sort_values('modeltype')
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
newdefdr['percent_activated2d'] = np.nan
# go through every row and for each fiber find the 2D activation percent
for row in newdefdr.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdefdr.query(
        'modeltype == "3DM" and nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and sim == @row.sim and master_fiber_index == @row.master_fiber_index and deformation==@row.deformation'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert val is not np.nan
    newdefdr.loc[row.Index, 'percent_activated3d'] = val
# %% plot act order
newdefdr['contact'] = newdefdr['contact'].replace(
    {'3D': "3DM", "cathodic": "cathodic-leading\n2DEM", "anodic": "anodic-leading\n2DEM", "center": "center\n2DEM"}
)
newdefdr['contact'] = pd.Categorical(
    newdefdr['contact'],
    categories=["3DM", "cathodic-leading\n2DEM", "center\n2DEM", "anodic-leading\n2DEM"],
    ordered=True,
)
# plot
sns.set(font_scale=1.25, style='whitegrid')
plt.figure()
g = sns.relplot(
    kind='scatter',
    row='contact',
    col='nerve_label',
    data=newdefdr.query("fiber_diam in [3]"),
    y='percent_activated3d',
    x='percent_activated',
    palette='plasma',
    hue='percent_activated3d',
    # hue='inner',
    facet_kws={'margin_titles': True},
    s=25,
)
plt.gca().set_aspect('equal')
# %% dose-response
panelA = defdr.query(f"nerve_label in {defsamples} and contact in {main_comparison} and deformation=='Undeformed'")
panelA['panel'] = 'Undeformed 2DEM\n(center slice) vs. 3DM'

panelB = defdr.query(f"modeltype=='3DM' and deformation!='3D-3D' and contact in {main_comparison}")
panelB['panel'] = '3DM Undeformed\nvs. Deformed'

panelC = defdr.query(f"deformation!='Undeformed' and contact in {main_comparison}")
panelC['panel'] = 'Deformed 2DEM\n(center slice) vs. 3DM'

panelD = defdr.query(f"deformation!='Undeformed' and contact in {cath_comparison}")
panelD['panel'] = 'Deformed 2DEM (cathodic-leading\nslice) vs. 3DM'

alldr = pd.concat([panelA, panelB, panelC, panelD]).query('nerve_label=="2L"')
alldr['nerve_label'] = pd.Categorical(alldr['nerve_label'], categories=['2L'])

alldr['deformation'] = pd.Categorical(alldr['deformation'], categories=['Undeformed', '2D-3D', '3D-3D'], ordered=True)

sns.set(font_scale=1.75, style='whitegrid')

plt.figure()
g = sns.relplot(
    kind='line',
    style='deformation',
    data=alldr.query("fiber_diam in [3]"),
    y='percent_activated',
    x='threshold',
    units='nerve_label',
    hue='modeltype',
    palette=pal2d3d,
    estimator=None,
    linewidth=4,
    facet_kws={'sharex': True, 'margin_titles': False},
    row='panel',
)

g.legend.set_title('')

for ax in g.axes.ravel():
    ax.set_ylim([0, 1])
    ax.set_xlim([0, None])
rownames(g, row_template="Proportion Activated\n")
g.set_titles(row_template='{row_name}', col_template="")
# plt.subplots_adjust(hspace=0.25)

g.set_xlabels('Threshold (mA)')
# change the line width for the legend
for line, l in zip(g.legend.get_lines(), g.legend.get_texts()):
    line.set_linewidth(4.0)
    if l.get_text() in ['deformation', 'modeltype']:
        l.set_text('')
g.fig.set_size_inches(8, 16)

# %%
sns.displot(
    kind='hist',
    x='threshold',
    col='nerve_label',
    row='contact',
    hue='deformation',
    data=defdr.query('fiber_diam==13'),
    facet_kws={'sharex': True, 'margin_titles': True},
)
# plt.gca().set_xscale('log')

import statsmodels.api as sm

# %%
from statsmodels.formula.api import ols

model = ols('threshold ~ nerve_label', data=defdr).fit()
print(model.summary())
anova_table = sm.stats.anova_lm(model, typ=2)
print(anova_table)
defdr.to_csv(r'C:\Users\dpm42\My Drive\ClassFilesGrad\SLOTKIN\analysis\defdr.csv', index=False)
matched.to_csv(r'C:\Users\dpm42\My Drive\ClassFilesGrad\SLOTKIN\analysis\matched.csv', index=False)

# datafromr= pd.read_csv(r'C:\Users\dpm42\My Drive\ClassFilesGrad\SLOTKIN\analysis\quantiled_dr.csv')
# out = AnovaRM(data=datafromr, depvar='log_quantile_threshold',
#               subject='nerve_label', within=['contact','deformation','fiber_diam','modeltype','quantile'],
#               typ=1).fit()
# print(out)
# %% LCCC comparison maincomp


def concordance_correlation_coefficient(y_true, y_pred):
    """Concordance correlation coefficient."""
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    # Raw data
    dct = {'y_true': y_true, 'y_pred': y_pred}
    df = pd.DataFrame(dct)
    # Remove NaNs
    df = df.dropna()
    # Pearson product-moment correlation coefficients
    y_true = df['y_true']
    y_pred = df['y_pred']
    cor = np.corrcoef(y_true, y_pred)[0][1]
    # Means
    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)
    # Population variances
    var_true = np.var(y_true)
    var_pred = np.var(y_pred)
    # Population standard deviations
    sd_true = np.std(y_true)
    sd_pred = np.std(y_pred)
    # Calculate CCC
    numerator = 2 * cor * sd_true * sd_pred
    denominator = var_true + var_pred + (mean_true - mean_pred) ** 2

    return numerator / denominator


def CCC_wp(y_true, y_pred):
    """Concordance correlation coefficient.

    With pval
    """
    r = concordance_correlation_coefficient(y_true, y_pred)
    t = r * np.sqrt(len(y_true) - 2) / np.sqrt(1 - r**2)
    p = 2 * (1 - stats.t.cdf(abs(t), len(y_true) - 2))
    return r, p


nsimdata = concats.query(
    f'fiber_diam in [3,13] and contact in {main_comparison}'
)  # TODO replace all cath comparison with non
g = sns.relplot(
    data=nsimdata.rename(columns={'nerve_label': 'Sample'}),
    kind='scatter',
    col='fiber_diam',
    x='threshold',
    y='threshold3d',
    # hue='Sample',
    color='white',
    row='deformation',
    # s=20,
    palette='colorblind',
    facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
    edgecolor='black',
    linewidth=1,
    alpha=1,
)
# TODO Clean up this calc
for diam, pos, ax, deformation in zip(
    [3, 13, 3, 13], (0.2, 0.8, 0.2, 0.8), g.axes.ravel(), ["Undeformed", "Undeformed", "3D-3D", "3D-3D"]
):
    rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
    r, p = CCC_wp(rdata.threshold3d, rdata.threshold)
    perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    lim = ax.get_xlim()[1]
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    print(f'{diam} {deformation} μm: {r ** 2:.2f} - p={p}')
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim], [0, lim], '--k', linewidth=2, label='1:1 line')
    # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])

    # ax.set_yticks(ax.get_xticks())
g.set_titles(col_template='D: {col_name} μm')
g.set_xlabels('2DEM Threshold (mA)')
g.set_ylabels('3DM Threshold (mA)')
# %% activation order
newdefdr = defdr.copy().query('deformation!="2D-3D"').reset_index()
newdefdr = newdefdr.sort_values('modeltype')
sns.set(font_scale=1.25)
sns.set_style('whitegrid')
newdefdr['percent_activated2d'] = np.nan
# go through every row and for each fiber find the 2D activation percent
for row in newdefdr.itertuples():
    # find the 2D threshold for this fiber (same nerve, fiber diameter, and master fiber index)
    thisdat = newdefdr.query(
        'modeltype == "3DM" and nerve_label == @row.nerve_label and fiber_diam == @row.fiber_diam and sim == @row.sim and master_fiber_index == @row.master_fiber_index and deformation==@row.deformation'
    )
    assert len(thisdat) == 1
    val = thisdat.percent_activated.values[0]
    assert val is not np.nan
    if row.modeltype == "3DM":
        assert newdefdr.loc[row.Index, 'percent_activated'] == val
    newdefdr.loc[row.Index, 'percent_activated3d'] = val
newdefdr['contact'] = newdefdr['contact'].replace(
    {'3D': "3DM", "cathodic": "cathodic-leading\n2DEM", "anodic": "anodic-leading\n2DEM", "center": "center\n2DEM"}
)
newdefdr['contact'] = pd.Categorical(
    newdefdr['contact'],
    categories=["3DM", "cathodic-leading\n2DEM", "center\n2DEM", "anodic-leading\n2DEM"],
    ordered=True,
)
# %% plot
sns.set(font_scale=1.75, style='whitegrid')
plt.figure()
g = sns.relplot(
    kind='scatter',
    data=newdefdr.query("nerve_label=='5R' and fiber_diam==13 and deformation=='3D-3D'"),
    y='percent_activated',
    x='percent_activated3d',
    palette='colorblind',
    hue='contact',
    facet_kws={'margin_titles': True},
    s=25,
)
for ax in g.axes.ravel():
    ax.set_aspect('equal')
    lim = max([ax.get_ylim()[1], ax.get_xlim()[1]])
    ax.plot([0, lim], [0, lim], 'k--')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
g.set_xlabels('Prop. 3DM fibers\nactivated')
g.set_ylabels('Prop. fibers activated')
g.set_titles(col_template='{col_name}')
newdefdr.to_csv(r'C:\Users\dpm42\My Drive\ClassFilesGrad\SLOTKIN\analysis\actorder.csv', index=False)
# %%
sns.set(font_scale=2, style='whitegrid')
datproc = (
    concats.query("deformation!='2D-3D'")
    .groupby(["deformation", "fiber_diam", "contact", "nerve_label", "inner"])[["threshold", "threshold3d"]]
    .agg(np.mean)
    .dropna()
    .reset_index()
)
g = sns.relplot(
    data=datproc.query('fiber_diam in [3,13]'),
    kind="scatter",
    row="fiber_diam",
    col="contact",
    x="threshold",
    y="threshold3d",
    hue='nerve_label',
    s=100,
    style='deformation',
    facet_kws=dict(sharex=False, sharey=False, margin_titles=True),
)
for ax in g.axes.ravel():
    limmax = max([ax.get_ylim()[1], ax.get_xlim()[1]])
    ax.set_ylim([0, limmax])
    ax.set_xlim([0, limmax])
    ax.plot([0, limmax], [0, limmax], 'k--', linewidth=2, alpha=0.5)
g.set_ylabels('3DM Thresholld (mA)')
g.set_xlabels('2DEM Threshold (mA)')
g.set_titles(col_template='{col_name}', row_template='{row_name}')
g.legend.get_texts()[0].set_text('Sample')
# %% new r2
sns.reset_orig()
sns.set(style='whitegrid', font_scale=1.5)
nsimdata = concats.query('fiber_diam in [3,13]')  # TODO replace all cath comparison with non
nsimdata['contact'] = pd.Categorical(nsimdata['contact'], categories=['anodic', 'center', 'cathodic'], ordered=True)
nsimdata['slice'] = nsimdata['contact'].copy()
nsimdata['deformation'].replace({'3D-3D': 'Structural'}, inplace=True)
g = sns.relplot(
    data=nsimdata.rename(columns={'nerve_label': 'Sample'}),
    kind='scatter',
    row='fiber_diam',
    x='threshold',
    y='threshold3d',
    # hue='Sample',
    color='white',
    hue='deformation',
    col='slice',
    # s=20,
    palette=['#00BFC4', '#F8766D'],
    facet_kws={'sharex': False, 'sharey': False, 'margin_titles': True},
    edgecolor='black',
    linewidth=1,
    alpha=1,
)
g.set_titles(col_template='{col_name}', row_template='{row_name}')


g.set_xlabels('2DEM Threshold (mA)')
g.set_ylabels('3DM Threshold (mA)')
for ax in g.axes.ravel():
    # rdata = nsimdata.query(f'fiber_diam=={diam} and deformation=="{deformation}"')
    # r,p = CCC_wp(rdata.threshold3d, rdata.threshold)
    # perc = sum(rdata.threshold > rdata.threshold3d) / len(rdata.threshold)
    lim = ax.get_xlim()[1]
    # add correlation to plot
    # ax.text(0.02, 0.91, f'$R^2={r**2:.2f}$', transform=ax.transAxes)
    # print(f'{diam} {deformation} μm: {r ** 2:.2f} - p={p}')
    # ax.set_title(f'{diam} μm')
    ax.plot([0, lim], [0, lim], '--k', linewidth=2, label='1:1 line')
    # ax.set_aspect('equal', 'box')
    # ax.apply_aspect()
    ax.set_xlim([0, lim])
    ax.set_ylim([0, lim])
