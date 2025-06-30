def plot_kmer_hist(kmer_collect, save_path=None, top_n=20, kmer_list=None, species_colors=None):
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    # Get species in desired order
    if species_colors:
        species = tuple(species_colors.keys())
    else:
        species = tuple(kmer_collect.kmer_mat_idc_df["species"].unique().to_list())

    # Compute mean kmer counts per species
    kmer_means = {}
    for sp in species:
        n_s = kmer_collect.kmer_mat_idc_df.filter(pl.col("species") == sp)["start index"].min()
        n_e = kmer_collect.kmer_mat_idc_df.filter(pl.col("species") == sp)["end index"].max()
        print(f"{sp} indices: {n_s}-{n_e}")
        kmer_count = kmer_collect[n_s:n_e].kmer_hist
        kmer_means[sp] = kmer_count.mean(axis=0)

    # Get kmer names and all mean values
    kmer_names = kmer_collect.get_kmer_column_names()
    species_vals = [kmer_means[sp] for sp in species]
    stacked = np.vstack(species_vals).T  # shape: (num_kmers, num_species)

    # Determine which kmers to plot
    if kmer_list is not None:
        kmer_name_to_index = {name: i for i, name in enumerate(kmer_names)}
        missing = [k for k in kmer_list if k not in kmer_name_to_index]
        if missing:
            raise ValueError(f"The following kmers are not present: {missing}")
        kmer_indices = [kmer_name_to_index[k] for k in kmer_list]
        selected_kmer_names = kmer_list
        selected_kmer_vals = stacked[kmer_indices]
    else:
        total_counts = stacked.sum(axis=1)
        top_indices = np.argsort(total_counts)[::-1][:top_n]
        selected_kmer_names = [kmer_names[i] for i in top_indices]
        selected_kmer_vals = stacked[top_indices]

    # Plot setup
    num_kmers = len(selected_kmer_names)
    x = np.arange(num_kmers)
    width = 0.8 / len(species)
    fig, ax = plt.subplots(figsize=(max(10, num_kmers), 6))

    for i, sp in enumerate(species):
        vals = selected_kmer_vals[:, i]
        offset = i * width
        color = species_colors.get(sp) if species_colors else None
        rects = ax.bar(x + offset, vals, width, label=sp, color=color)

    ax.set_xticks(x + width * (len(species) - 1) / 2)
    ax.set_xticklabels(selected_kmer_names, rotation=90, ha="center", fontsize=9)

    ax.set_xlabel("K-mer")
    ax.set_ylabel("Mean k-mer count")
    ax.set_title(
        f"{'Selected' if kmer_list else f'Top {top_n}'} k-mers by mean count (grouped by k-mer, colored by species)"
    )
    ax.legend(title="Species", bbox_to_anchor=(1.01, 1), loc="upper left")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches="tight")
    plt.show()

    return selected_kmer_vals, selected_kmer_names
