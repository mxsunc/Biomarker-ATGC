import pylab as plt
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

D, samples = pickle.load(open(cwd / 'sim_data' / 'regression' / 'experiment_1' / 'sim_data.pkl', 'rb'))

idx_test, instance_sum_predictions = pickle.load(open(cwd / 'sim_data' / 'regression' / 'experiment_1' / 'instance_model_sum_predictions.pkl', 'rb'))
idx_test, instance_mean_predictions = pickle.load(open(cwd / 'sim_data' / 'regression' / 'experiment_1' / 'instance_model_mean_predictions.pkl', 'rb'))
idx_test, sample_sum_predictions = pickle.load(open(cwd / 'sim_data' / 'regression' / 'experiment_1' / 'sample_model_sum_predictions.pkl', 'rb'))
idx_test, sample_mean_predictions = pickle.load(open(cwd / 'sim_data' / 'regression' / 'experiment_1' / 'sample_model_mean_predictions.pkl', 'rb'))

##get x_true and y_true

x_true = []
for sample_idx in idx_test:
    variants = D['class'][np.where(D['sample_idx'] == sample_idx)]
    x_true.append(len(np.where(variants != 0)[0]))

y_true = np.array(samples['values'])[idx_test]

predictions = [y_true, instance_sum_predictions, instance_mean_predictions, sample_sum_predictions, sample_mean_predictions]


z_order = np.arange(0, len(x_true))
colors = ['k', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(
top=1.0,
bottom=0.075,
left=0.035,
right=1.0,
hspace=0.2,
wspace=0.2)
for index, i in enumerate(['True Value', 'Instance Sum', 'Instance Mean', 'Sample Sum', 'Sample Mean']):
    if index != 0:
        for index2, (x, y) in enumerate(zip(x_true, np.exp(predictions[index][0]) - 1)):
            if index2 == len(x_true) - 1:
                ax.scatter(x, y, color=colors[index], edgecolor='k', zorder=np.random.choice(z_order), label=i)
            else:
                ax.scatter(x, y, color=colors[index], edgecolor='k', zorder=np.random.choice(z_order))

    else:
        for index2, (x, y) in enumerate(zip(x_true, predictions[index])):
            if index2 == len(x_true) - 1:
                ax.scatter(x, y, color=colors[index], zorder=np.random.choice(z_order), label=i)
            else:
                ax.scatter(x, y, color=colors[index], zorder=np.random.choice(z_order))

ax.set_yticks([])
ax.set_xticks([])

ax.set_ylabel('Bag Value', fontsize=24, labelpad=-10)
ax.set_xlabel('Key Instance Count', fontsize=24, labelpad=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.legend(frameon=False, loc='upper left', fontsize=14)

plt.savefig(cwd / 'sim_data' / 'regression' / 'experiment_1' / 'figure_predictions.png', dpi=300)

