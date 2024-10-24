
import os
import sys
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# read in data
data = pd.read_csv("Gantt_12.csv", parse_dates=['Start', 'End'], dayfirst=True)

# convert data to data frame
gant_df = pd.DataFrame(data)

# convert date strings into datetime objects
gant_df['Start'] = pd.to_datetime(gant_df['Start'])
gant_df['End'] = pd.to_datetime(gant_df['End'])

# sort dataFrame according to start date
gant_df = gant_df.sort_values(by='Start', ascending=False)

# plotting the Gant chart
fig, ax = plt.subplots(figsize=(8, len(gant_df) * 0.5))

# define colors for each of the tasks

custom_colors = [
    '#1f77b4',  # blue
    '#ff7f0e',  # orange
    '#2ca02c',  # green
    '#d62728',  # red
    '#9467bd',  # purple
    '#8c564b',  # brown
    '#e377c2',  # pink
    '#7f7f7f',  # gray
    '#bcbd22',  # olive
    '#17becf',  # cyan
    '#aec7e8',  # light blue
    '#ffbb78',  # light orange
    '#98df8a',  # light green
    '#ff9896',  # light red
    '#c5b0d5',  # light purple
    '#c49c94',  # light brown
]

for i, row in gant_df.iterrows():
    ax.barh(row['Task'], row['End'] - row['Start'], left=row['Start'],height=0.8, align='center', alpha=0.7)

# improve plot
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
ax.yaxis.grid(True)
ax.xaxis.grid(True)
ax.set_ylabel('Tasks')
ax.set_xlabel('Timeline_Dates')
ax.set_title('12 Months Gantt Chart')

plt.yticks(range(len(gant_df)), gant_df['Task'], fontsize=9)
plt.xticks(rotation=45)
plt.grid(True, which='both', axis='x', alpha=0.15)

plt.show()



