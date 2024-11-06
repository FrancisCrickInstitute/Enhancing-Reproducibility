import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

plt.rcParams['font.size'] = 26
plt.rcParams['axes.linewidth'] = 2

# Set the seed for reproducibility
np.random.seed(5)

# Define the threshold value for x
threshold = 2000
start = 2000

# Define the mean and standard deviation for the Gaussian distribution
# mean1 = -0.05
# mean1 = 1010
# mean2 = 990
mean1 = -0.02
mean2 = 0.1
mean3 = -0.1
# mean2 = 0.05
std_dev = 1.8

# Initialize the samples as empty arrays
sample1 = np.random.normal(mean1, std_dev, start)
sample2 = np.random.normal(mean2, std_dev, start)
sample3 = np.random.normal(mean3, std_dev, start)

# Initialize lists to store the sample sizes and p-values
sample_sizes = []
p_values = []

# Loop over x from 1 to the threshold value
for x in range(start, threshold + 1):
    # Generate new samples of 1 datapoint from the Gaussian distribution
    new_sample1 = np.random.normal(mean1, std_dev, 1)
    new_sample2 = np.random.normal(mean2, std_dev, 1)
    new_sample3 = np.random.normal(mean3, std_dev, 1)

    # Add the new samples to the existing samples
    sample1 = np.append(sample1, new_sample1)
    sample2 = np.append(sample2, new_sample2)
    sample3 = np.append(sample3, new_sample3)

    # Perform a t-test on the two samples
    t_stat1, p_val1 = stats.ttest_ind(sample1, sample2)
    t_stat2, p_val2 = stats.ttest_ind(sample1, sample3)
    t_stat3, p_val3 = stats.ttest_ind(sample2, sample3)
    print(f'P-Value 1: {p_val1}')
    print(f'P-Value 2: {p_val2}')
    print(f'P-Value 3: {p_val3}')

    # Combine the two samples into a single array
    data = np.concatenate([sample1, sample2, sample3])

    # Create an array of labels for the two samples
    labels = np.concatenate([np.full(len(sample1), 'Slide 1'), np.full(len(sample2), 'Slide 2'),
                             np.full(len(sample3), 'Slide 3')])

    # Create a DataFrame with the data and labels
    df = pd.DataFrame({'Data': data, 'Label': labels})

    plt.figure(figsize=(14, 10))

    # Create a swarmplot of the two distributions, with a boxplot overlaid
    sns.boxplot(x='Label', y='Data', data=df, boxprops=dict(facecolor='none', zorder=2),
                whiskerprops=dict(color="black", linewidth=2, zorder=2),
                capprops=dict(color="black", linewidth=2, zorder=2),
                medianprops=dict(color="black", linewidth=2, zorder=2),
                showfliers=False)
    sns.swarmplot(x='Label', y='Data', data=df, hue='Label', palette=['blue', 'gray', 'orange'], size=3, alpha=0.9,
                  zorder=1)

    # Annotate the plot with the result of the t-test
    # title_color = 'black' if p_val >= 0.05 else 'red'
    # plt.title(f'p-value: {p_val:.4f}', color=title_color)
    plt.xlabel('')
    plt.ylabel('Metric')
    #plt.ylabel(r'Mean Cell Area ($\mu$m$^2$)')
    #plt.ylim(bottom=600, top=1400)
    plt.tight_layout()

    # Save the plot to disk
    plt.savefig(f'./outputs/plots/anim/t_test_plot_{x}.png')

    # Clear the plot for the next iteration
    plt.close()

    # Append the sample size and p-value to the lists
    sample_sizes.append(x)
    #p_values.append(p_val)

    # plt.figure(figsize=(14, 10))
    #
    # # Create a line plot of p-value versus sample size
    # plt.scatter(sample_sizes, p_values, s=100)
    # plt.axhline(y=0.05, color='red', linestyle='dotted', label='p = 0.05')
    # plt.xlabel('Number of Cells')
    # plt.ylabel('P-value')
    # plt.xlim(left=start, right=threshold)
    # plt.ylim(top=1.0, bottom=0.0)
    # plt.tight_layout()
    # plt.savefig(f'./outputs/plots/anim2/p_value_vs_sample_size_{x}.png')
    #
    # # Clear the plot for the next iteration
    # plt.close()
