
# Add warning of number of cohorts per gene
df_counts = df_drivers.groupby("SYMBOL", as_index=False).agg({"COHORT": "count"})
df_counts.rename(columns={"COHORT": "num_cohorts"}, inplace=True)
df_drivers = df_drivers.merge(df_counts)
df_drivers["Warning_num_cohorts"] = df_drivers.apply(lambda row: True if row["num_cohorts"] == 1 else False, axis=1)


# TODO role and clusters 2D aminoacids


df = df[(df['P'] < PVALUE_THRESHOLD)]
if df.shape[0] > 0:
    # Get the amino acid coordinates
    df[["AA_START", "AA_END"]] = df.apply(
        lambda row: get_position_aa(reader, row["COORDINATES"], row["CHROMOSOME"], row["ENSID"]), axis=1)
    df["2D_CLUSTERS"] = df["AA_START"] + ":" + df["AA_END"]
    df = df.groupby(["SYMBOL", "COHORT"], as_index=False).agg(
        {"2D_CLUSTERS": lambda x: ','.join(set(map(str, x)))})
    clusters.append(df.copy())
