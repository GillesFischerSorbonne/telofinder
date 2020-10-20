from telofinder import analyze_telom_length as atl
import matplotlib.pyplot as plt

df, telo_df, telo_df_merged  = atl.run_telofinder("../telofinder/data/AFH_chrI.fasta", 0.8, 0.8, 8000)
df_res = df.reset_index()
df_W = df_res['level_3'] == "W"

ax1 = (
    df_res[df_W].loc[:, ["polynuc", "entropy"]]
    .plot()
    .legend(loc="center left", bbox_to_anchor=(0.75, 0.15))
)
ax1.set_title("AFH chr1 left")
plt.ylim(top=1.6)
plt.hlines(0.8, xmin = 0, xmax = 8000, colors='green', linestyles='dotted')
plt.text(3000, 0.75, 'Threshold', ha ='left', va ='top', color='green')
plt.hlines(1.45, xmin = 0, xmax = 246, colors='green', linestyles='solid')
plt.text(0, 1.55, 'Terminal Telomere', ha ='left', va ='top', color='green')
plt.hlines(1.45, xmin = 7080, xmax = 7160, colors='green', linestyles='solid')
plt.text(6000, 1.55, 'Internal Telomere', ha ='left', va ='top', color='green')
plt.title("Telomere detection in the first 8kb", loc="left")

ax2 = (
    df_res[df_W].loc[:, ["polynuc", "entropy"]]
    .plot()
    .legend(loc="center left", bbox_to_anchor=(0.75, 0.15))
)
ax2.set_title("AFH chr1 left")
plt.ylim(top=1.6)
plt.xlim(0, 600)
plt.hlines(0.8, xmin = 0, xmax = 600, colors='green', linestyles='dotted')
plt.text(300, 0.75, 'Threshold', ha ='left', va ='top', color='green')
plt.hlines(1.45, xmin = 0, xmax = 226, colors='green', linestyles='solid')
plt.text(0, 1.55, 'Terminal Telomere', ha ='left', va ='top', color='green')
plt.title("Zoom in terminal region [0-600bp]", loc="left")

ax3 = (
    df_res[df_W].loc[:, ["polynuc", "entropy"]]
    .plot()
    .legend(loc="center left", bbox_to_anchor=(0.75, 0.15))
)
ax3.set_title("AFH chr1 left")
plt.ylim(top=1.6)
plt.xlim(7000, 7300)
plt.hlines(0.8, xmin = 7000, xmax = 7300, colors='green', linestyles='dotted')
plt.text(7025, 0.75, 'Threshold', ha ='left', va ='top', color='green')
plt.hlines(1.45, xmin = 7089, xmax = 7135, colors='green', linestyles='solid')
plt.text(7089, 1.55, 'Internal Telomere', ha ='left', va ='top', color='green')
plt.title("Zoom in internal region [7000-7300bp]", loc="left")