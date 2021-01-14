import pylab as plt
import seaborn as sns
import pandas as pd
import numpy as np
import pickle
import pathlib
path = pathlib.Path.cwd()

if path.stem == 'ATGC2':
    cwd = path
else:
    cwd = list(path.parents)[::-1][path.parts.index('ATGC2')]
    import sys
    sys.path.append(str(cwd))

D, samples = pickle.load(open(cwd / 'sim_data' / 'classification' / 'experiment_1' / 'sim_data.pkl', 'rb'))

idx_test, mean_attentions = pickle.load(open(cwd / 'sim_data' / 'classification' / 'experiment_1' / 'mean_attentions.pkl', 'rb'))
idx_test, sum_attentions = pickle.load(open(cwd / 'sim_data' / 'classification' / 'experiment_1' / 'sum_attentions.pkl', 'rb'))
idx_test, both_attentions = pickle.load(open(cwd / 'sim_data' / 'classification' / 'experiment_1' / 'both_attentions.pkl', 'rb'))
idx_test, dynamic_attentions = pickle.load(open(cwd / 'sim_data' / 'classification' / 'experiment_1' / 'dynamic_attentions.pkl', 'rb'))


idx_test = idx_test[:20]
indexes = [np.where(D['sample_idx'] == idx) for idx in idx_test]

classes = []
for i in indexes:
    classes.append(D['class'][i])

types = np.concatenate(classes).shape[0] * [0] + np.concatenate(classes).shape[0] * [1] + np.concatenate(classes).shape[0] * [2] + np.concatenate(classes).shape[0] * [3]
classes = np.concatenate([np.concatenate(classes), np.concatenate(classes), np.concatenate(classes), np.concatenate(classes)])
attention = np.concatenate([np.concatenate(mean_attentions[2][:20]),
                            np.concatenate(sum_attentions[2][:20]),
                            np.concatenate(both_attentions[2][:20]),
                            np.concatenate(dynamic_attentions[2][:20])])

instance_df = pd.DataFrame({'attention': attention.flat, 'class': classes, 'type': types})

fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(top=0.927,
bottom=0.104,
left=0.043,
right=1.0,
hspace=0.2,
wspace=0.2)
sns.stripplot(x="type", y="attention", hue='class', data=instance_df, edgecolor='k', linewidth=1, jitter=.4)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_yticks([])
ax.tick_params(length=0, pad=-5)
ax.set_xticklabels(['Mean', 'Sum', 'Both', 'Dynamic'], fontsize=14)
ax.set_ylabel('Attention', fontsize=16)
ax.set_xlabel('Model', fontsize=16, labelpad=10)
ax.legend(frameon=False, loc=[.05, .9], title='Instance Class')
plt.savefig(cwd / 'sim_data' / 'classification' / 'experiment_1' / 'attention.png', dpi=300)

