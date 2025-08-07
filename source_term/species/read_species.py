import pandas as pd

# filename = 'run_dep_out_species0.csv'
filename = 'sub_species_out.csv'
df = pd.read_csv(filename, index_col=False)
isotope = 'Xe127'
isotope_keys = [key for key in df.keys() if isotope in key]
subdf = df[isotope_keys]
# subdf.to_csv('out.csv', index=False)
print(df[f'{isotope}_integral_no_removal'])
print(df[f'{isotope}_initial_integral'])
print(df[f'{isotope}_integral'])
# print(df[f'{isotope}_post'])
print(df[f'{isotope}_gas'])
print(df[f'{isotope}_fission_production'])
print(df[f'{isotope}_th_removal'])
print(df[f'{isotope}_mass_transfer_max'])
print(df['Xe_removal'])
print(df['fiss_rate_density'])