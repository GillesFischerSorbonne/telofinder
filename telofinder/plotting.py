

def plot_telom(telom_df):
    """Plotting the telomere detection on both left and right chromosome ends
    """
    df = telom_df.reset_index()
    for strand in ["W", "C"]:
        ax = (
            df.query("level_3==@strand")
            .loc[:, ["polynuc", "entropy", "predict_telom"]]
            .plot()
            .legend(loc="center left", bbox_to_anchor=(1, 0.5))
        )
        ax.set_title(strand)



