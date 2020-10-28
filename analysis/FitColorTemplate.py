# ------------------------------------------------------------------------------------------------------------------- #
# Read Other Type II SNe Data From The SN Archive Folder
# Read B-V Color Template (Stritzinger et al. 2018)
# ------------------------------------------------------------------------------------------------------------------- #
list_files = group_similar_files('', DIR_SNe + 'LC_Data/*.asc', exceptions='SWIFT')

colbv_df = pd.read_csv(os.path.join(DIR_TEMP, file_bvcolor), sep='\s+', names=cols_temp)
colvr_df = pd.read_csv(os.path.join(DIR_TEMP, file_vrcolor), sep='\s+', names=cols_temp)
colri_df = pd.read_csv(os.path.join(DIR_TEMP, file_ricolor), sep='\s+', names=cols_temp)

display_text("Loaded B-V, V-r and r-i Templates from Stritzinger et al. (2018)")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Copy Pandas DataFrame Into Subset DataFrames Required For Various Plots
# ------------------------------------------------------------------------------------------------------------------- #
rawopt_df = pd.read_csv(DIR_SNe + 'LC_Data/' + file_opt, sep='\s+', comment='#')
rawopt_df = rawopt_df.replace('INDEF', np.nan).astype('float64')
rawopt_df['Phase'] = rawopt_df['JD'] - date_bmax

# outputopt_df = multicol_to_fluxdf(name_SN, rawopt_df)
# vabs_df = outputopt_df[outputopt_df['FILTER'] == 'V'].copy()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Optical Colours
# ------------------------------------------------------------------------------------------------------------------- #
dataopt_df = multicol_to_fluxdf(name_SN, rawopt_df).set_index('Date')

snmag_df = calc_extcorcolorframe(unorgframe_to_orgframe(dataopt_df, column='FMAG'))
snerr_df = calc_extcorcolorframe(unorgframe_to_orgframe(dataopt_df, column='FERR'), err=True)
snGalmag_df = calc_extcorcolorframe(unorgframe_to_orgframe(dataopt_df, column='FMAG'), ebv=EBVMW_mag)
snGalerr_df = calc_extcorcolorframe(unorgframe_to_orgframe(dataopt_df, column='FERR'), ebverr=EBVMW_err, err=True)
# tmag_df = calc_extcorcolorframe(unorgframe_to_orgframe(dataopt_df, column='FMAG'), ebv=0.26)
# terr_df = calc_extcorcolorframe(unorgframe_to_orgframe(dataopt_df, column='FERR'), ebverr=0.05, err=True)

max_epoch = 2458200.0
snmag_df = snmag_df[snmag_df['Phase'] < max_epoch].set_index('Phase')
snerr_df = snerr_df[snerr_df['Phase'] < max_epoch].set_index('Phase')
snGalmag_df = snGalmag_df[snGalmag_df['Phase'] < max_epoch].set_index('Phase')
snGalerr_df = snGalerr_df[snGalerr_df['Phase'] < max_epoch].set_index('Phase')
# tmag_df = tmag_df[tmag_df['Phase'] < max_epoch].set_index('Phase')
# terr_df = terr_df[terr_df['Phase'] < max_epoch].set_index('Phase')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Power Law Fit To The Pre-Maximum Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

for file_name in list_files:
    name = file_name.split('/')[-1].split('.')[0]
    if name in ['iPTF13bvn', '2012au', '2009jf', '2007Y']:
        datacomp_df = pd.read_csv(file_name, sep='\s+', comment='#')
        if 'Phase' not in datacomp_df.columns:
            datacomp_df['Phase'] = datacomp_df['JD'] - sndata_df.loc[name, 'DateBMax']
        if 'Date' not in datacomp_df.columns:
            datacomp_df['Date'] = datacomp_df['JD'].apply(lambda x: jd_to_cald(x))
        datacomp_df = datacomp_df.replace('INDEF', np.nan).sort_values(by='Phase')
        plot_color(ax, datacomp_df, name, 'B-V')

# ax.errorbar(col_df['Phase'], col_df['Color'], yerr=col_df['ColorErr'], marker='', c='dimgrey', ls='',
#             capsize=3, capthick=0.5, label='_nolegend_')
ax.plot(colbv_df['Phase'], colbv_df['Color'], marker='', c='k', ls='-', lw=2, ms=15, label='B-V Template (S18)')
ax.plot(colbv_df['Phase'], colbv_df['Color'] - colbv_df['ColorErr'], marker='', c='k', ls='--', lw=1,
        label='_nolegend_')
ax.plot(colbv_df['Phase'], colbv_df['Color'] + colbv_df['ColorErr'], marker='', c='k', ls='--', lw=1,
        label='_nolegend_')
ax.fill_between(colbv_df['Phase'], colbv_df['Color'] - colbv_df['ColorErr'], colbv_df['Color'] + colbv_df['ColorErr'],
                color='dodgerblue')

plot_sncolor(ax, snGalmag_df, snGalerr_df, color='B-V', ms=12)
plot_sncolor(ax, snmag_df, snerr_df, color='B-V', c='r', marker='s', ms=12)
# plot_sncolor(ax, tmag_df, terr_df, color='B-V', c='g', marker='s', ms=12)

ax.set_xlim(-15, 25)
# ax.set_ylim(42.3, 42.7)
ax.legend(markerscale=1.4, fontsize=15)
set_plotparams(ax, yticks=(0.2, 0.02), xticks=(10, 1), fs=16)
ax.set_ylabel('B-V [mag]', fontsize=18)
ax.set_xlabel('Time Since B-Band Maximum [Days]', fontsize=18)

fig.savefig('PLOT_StritzingerColorB-V.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Power Law Fit To The Pre-Maximum Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

for file_name in list_files:
    name = file_name.split('/')[-1].split('.')[0]
    if name in ['iPTF13bvn', '2012au', '2009jf', '2007Y']:
        datacomp_df = pd.read_csv(file_name, sep='\s+', comment='#')
        if 'Phase' not in datacomp_df.columns:
            datacomp_df['Phase'] = datacomp_df['JD'] - sndata_df.loc[name, 'DateBMax']
        if 'Date' not in datacomp_df.columns:
            datacomp_df['Date'] = datacomp_df['JD'].apply(lambda x: jd_to_cald(x))
        datacomp_df = datacomp_df.replace('INDEF', np.nan).sort_values(by='Phase')
        plot_color(ax, datacomp_df, name, 'V-R')

snGalmag_df.index = snGalmag_df.index - 1.8
snGalerr_df.index = snGalerr_df.index - 1.8
snmag_df.index = snmag_df.index - 1.8
snerr_df.index = snerr_df.index - 1.8

# ax.errorbar(col_df['Phase'], col_df['Color'], yerr=col_df['ColorErr'], marker='', c='dimgrey', ls='',
#             capsize=3, capthick=0.5, label='_nolegend_')
ax.plot(colvr_df['Phase'], colvr_df['Color'], marker='', c='k', ls='-', lw=2, ms=15, label='V-r Template (S18)')
ax.plot(colvr_df['Phase'], colvr_df['Color'] - colvr_df['ColorErr'], marker='', c='k', ls='--', lw=1,
        label='_nolegend_')
ax.plot(colvr_df['Phase'], colvr_df['Color'] + colvr_df['ColorErr'], marker='', c='k', ls='--', lw=1,
        label='_nolegend_')
ax.fill_between(colvr_df['Phase'], colvr_df['Color'] - colvr_df['ColorErr'], colvr_df['Color'] + colvr_df['ColorErr'],
                color='dodgerblue')

plot_sncolor(ax, snGalmag_df, snGalerr_df, color='V-R', ms=12)
plot_sncolor(ax, snmag_df, snerr_df, color='V-R', c='r', marker='s', ms=12)
# plot_sncolor(ax, tmag_df, terr_df, color='V-R', c='g', marker='s', ms=12)

ax.set_xlim(-15, 40)
ax.set_ylim(-0.1, 0.73)
ax.legend(markerscale=1.3, loc=2, fontsize=15)
set_plotparams(ax, yticks=(0.1, 0.01), xticks=(10, 1), fs=16)
ax.set_ylabel('V-R [mag]', fontsize=18)
ax.set_xlabel('Time Since V-Band Maximum [Days]', fontsize=18)

fig.savefig('PLOT_StritzingerColorV-R.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Power Law Fit To The Pre-Maximum Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

for file_name in list_files:
    name = file_name.split('/')[-1].split('.')[0]
    if name in ['iPTF13bvn', '2012au', '2009jf', '2007Y']:
        datacomp_df = pd.read_csv(file_name, sep='\s+', comment='#')
        if 'Phase' not in datacomp_df.columns:
            datacomp_df['Phase'] = datacomp_df['JD'] - sndata_df.loc[name, 'DateBMax']
        if 'Date' not in datacomp_df.columns:
            datacomp_df['Date'] = datacomp_df['JD'].apply(lambda x: jd_to_cald(x))
        datacomp_df = datacomp_df.replace('INDEF', np.nan).sort_values(by='Phase')
        plot_color(ax, datacomp_df, name, 'R-I')

snGalmag_df.index = snGalmag_df.index - 1.4
snGalerr_df.index = snGalerr_df.index - 1.4
snmag_df.index = snmag_df.index - 1.4
snerr_df.index = snerr_df.index - 1.4

# ax.errorbar(col_df['Phase'], col_df['Color'], yerr=col_df['ColorErr'], marker='', c='dimgrey', ls='',
#             capsize=3, capthick=0.5, label='_nolegend_')
ax.plot(colri_df['Phase'], colri_df['Color'], marker='', c='k', ls='-', lw=2, ms=15, label='r-i Template (S18)')
ax.plot(colri_df['Phase'], colri_df['Color'] - colri_df['ColorErr'], marker='', c='k', ls='--', lw=1,
        label='_nolegend_')
ax.plot(colri_df['Phase'], colri_df['Color'] + colri_df['ColorErr'], marker='', c='k', ls='--', lw=1,
        label='_nolegend_')
ax.fill_between(colri_df['Phase'], colri_df['Color'] - colri_df['ColorErr'], colri_df['Color'] + colri_df['ColorErr'],
                color='dodgerblue')

plot_sncolor(ax, snGalmag_df, snGalerr_df, color='R-I', ms=12)
plot_sncolor(ax, snmag_df, snerr_df, color='R-I', c='r', marker='s', ms=12)
# plot_sncolor(ax, tmag_df, terr_df, color='R-I', c='g', marker='s', ms=12)

ax.set_xlim(-15, 40)
ax.set_ylim(-0.21, 0.61)
ax.legend(markerscale=1.3, loc=2, fontsize=15)
set_plotparams(ax, yticks=(0.1, 0.01), xticks=(10, 1), fs=16)
ax.set_ylabel('R-I [mag]', fontsize=18)
ax.set_xlabel('Time Since R-Band Maximum [Days]', fontsize=18)

fig.savefig('PLOT_StritzingerColorR-I.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
