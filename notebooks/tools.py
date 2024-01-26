def normalize_well_format(well):
    match = re.match(r"([A-Za-z])([0-9]+)", well, re.I)
    if match:
        items = match.groups()
        return f"{items[0]}{int(items[1]):02d}"
    return well


def load_and_prepare_data(file_path, plate_number):
    df = pd.read_csv(file_path)
    df = df[df['Plate'] == plate_number]
    df['Control Type'] = df['Control Type'].fillna('Treated').replace('', 'Treated')
    df['Well'] = df['Well'].apply(normalize_well_format)
    return df


# Define a function to calculate the 95% confidence interval
def confidence_interval_95(data):
    ci = scipy.stats.t.interval(0.95, len(data) - 1, loc=np.mean(data), scale=scipy.stats.sem(data))
    return ci


# Define a function to calculate the interquartile range
def iqr(data):
    return np.subtract(*np.percentile(data, [75, 25]))