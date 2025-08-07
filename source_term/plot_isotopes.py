import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

N_pow_steps = len(pd.read_csv('power_history.csv')) - 1

# Read the CSV file
cases = ['th_removal_fractions', 'thermo_removal_fractions', 'sum_sets']
labels = ['Removal Fractions', 'Removal Fractions', 'Sum Set Values']
for index, case in enumerate(cases):
    df_list = []
    for power_index in range(N_pow_steps):
        csv_file = f'run_dep_{case}_{power_index}.csv'  # Change this to your actual file name
        df = pd.read_csv(csv_file)
        if power_index == 0:
            df_list.append(df)
        else:
            df_list.append(df[1:])

    # Now concat these data files
    data = pd.concat(df_list, ignore_index=True)

    # Set the time column as the x-axis
    time = data['time']

    # Create a figure and axis for plotting
    plt.figure(figsize=(12, 8))

    # Define a list of markers and line styles
    markers = ['o', 's', 'D', '^', 'v', 'x', '*', '+', 'P', 'H']
    linestyles = ['-', '--', '-.', ':']

    # Create a color map for random colors
    colors = plt.cm.rainbow(np.linspace(0, 1, len(data.columns) - 1))

    # Plot each sum set with different colors, markers, and line styles
    for i, column in enumerate(data.columns[1:]):  # Skip the first column (time)
        plt.plot(time, data[column], 
            label=column, 
            color=colors[i], 
            marker=markers[i % len(markers)], 
            linestyle=linestyles[i % len(linestyles)])

    # Customize the plot
    plt.xlabel('Time (s)')
    plt.ylabel(labels[index])
    plt.legend()
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.6))
    if case != 'thermo_removal_fractions':
        plt.yscale('log')
    plt.tight_layout()

    # Show the plot
    plt.savefig(f'{case}.png', dpi=500)