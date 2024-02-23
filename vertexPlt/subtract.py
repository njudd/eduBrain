# Taking the mean vertex value of the top and bottom 15% of SES & PGS and subtracting them.

# https://github.com/njudd/imagen/blob/master/sebastian/CT_import_script.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import nibabel.freesurfer.io as fsio
import nibabel.freesurfer.mghformat as mgh
import random
import seaborn as sns
from surfer import Brain




from scipy.interpolate import griddata
import nibabel as freesurfer
import statistics as stats


## DATA SECTION
# importing df
imagen = pd.read_csv("~/Projects/imaging/imagen/data/IMAGEN_Master_20190808.csv")
imagen = imagen[imagen.subs_for_finalstudy == 1]  # selecting relevant cols

# got the data from Jeshua dir
lh = np.load("/Users/nicjud/Projects/imaging/imagen/nick/analysis/py_vectorSubtraction/data/TP1_Surf_lh_sm10_h.npy")
rh = np.load("/Users/nicjud/Projects/imaging/imagen/nick/analysis/py_vectorSubtraction/data/TP1_Surf_rh_sm10_h.npy")
# reloading for noNAN
lh_noNAN = np.load("/Users/nicjud/Projects/imaging/imagen/nick/analysis/py_vectorSubtraction/data/TP1_Surf_lh_sm10_h.npy")
rh_noNAN = np.load("/Users/nicjud/Projects/imaging/imagen/nick/analysis/py_vectorSubtraction/data/TP1_Surf_rh_sm10_h.npy")

# np.sum(lh == 0.0)/551 # 13887 rh = 13916
# np.sum(lh[1,:] == 0.0) # 13887 rh = 13916

# is loop makes lh & rh values of 0.0 NANs, need for density histograms

# this is making a mask for the flat plots
lh_mask = lh[1, :] < 0.1  # I had it ==0.0 yet it kept a little dot in the middle
rh_mask = rh[1, :] < 0.1

for fix in range(0, lh_noNAN.shape[0]):
    lh_noNAN[fix, :][lh_noNAN[fix, :] == 0.0] = np.nan
    rh_noNAN[fix, :][rh_noNAN[fix, :] == 0.0] = np.nan


# checking that its in the right order
imagen = imagen.reset_index()
# first sub is the same
# imagen.lh_WhiteSurfArea_area_MRI1_total[0]
# lh.sum(axis=1)[0]
# last sub is the same
# imagen.lh_WhiteSurfArea_area_MRI1_total[550]
# lh.sum(axis=1)[550]

# making a logical above and bellow one SD for SES (in a norm distro this would be 16% of subs in each category)
lowSES = imagen[imagen.all_SES < stats.mean(imagen.all_SES) - stats.stdev(imagen.all_SES)]
highSES = imagen[imagen.all_SES > stats.mean(imagen.all_SES) + stats.stdev(imagen.all_SES)]

# indexing the vertices by the SES dfs
lh_lowSES = lh[lowSES.index]  # left hemi
lh_highSES = lh[highSES.index]
rh_lowSES = rh[lowSES.index]  # right hemi
rh_highSES = rh[highSES.index]
# same thing for values without zeros
lh_lowSES_noNAN = lh_noNAN[lowSES.index]  # left hemi
lh_highSES_noNAN = lh_noNAN[highSES.index]
rh_lowSES_noNAN = rh_noNAN[lowSES.index]  # right hemi
rh_highSES_noNAN = rh_noNAN[highSES.index]

# now getting the mean per vertex
lh_lowSES = lh_lowSES.sum(axis=0)
lh_highSES = lh_highSES.sum(axis=0)
diff_SES_lh = lh_highSES - lh_lowSES
rh_lowSES = rh_lowSES.sum(axis=0) # right hemi
rh_highSES = rh_highSES.sum(axis=0)
diff_SES_rh = rh_highSES - rh_lowSES

# same thing for values without zeros
lh_lowSES_noNAN = lh_lowSES_noNAN.sum(axis=0)
lh_highSES_noNAN = lh_highSES_noNAN.sum(axis=0)
diff_SES_lh_noNAN = lh_highSES_noNAN - lh_lowSES_noNAN
rh_lowSES_noNAN = rh_lowSES_noNAN.sum(axis=0)  # right hemi
rh_highSES_noNAN = rh_highSES_noNAN.sum(axis=0)
diff_SES_rh_noNAN = rh_highSES_noNAN - rh_lowSES_noNAN

# making a logical above and bellow one SD for PGS (in a norm distro this would be 16% of subs in each category)
lowPGS = imagen[imagen.PGS < stats.mean(imagen.PGS) - stats.stdev(imagen.PGS)]
highPGS = imagen[imagen.PGS > stats.mean(imagen.PGS) + stats.stdev(imagen.PGS)]

# indexing the vertices by the SES dfs
lh_lowPGS = lh[lowPGS.index]
lh_highPGS = lh[highPGS.index]
rh_lowPGS = rh[lowPGS.index]  # right hemi
rh_highPGS = rh[highPGS.index]
# same thing without zeros
lh_lowPGS_noNAN = lh_noNAN[lowPGS.index]
lh_highPGS_noNAN = lh_noNAN[highPGS.index]
rh_lowPGS_noNAN = rh_noNAN[lowPGS.index]  # right hemi
rh_highPGS_noNAN = rh_noNAN[highPGS.index]

# now getting the mean per vertex
lh_lowPGS = lh_lowPGS.sum(axis=0)
lh_highPGS = lh_highPGS.sum(axis=0)
diff_PGS_lh = lh_highPGS - lh_lowPGS
rh_lowPGS = rh_lowPGS.sum(axis=0)  # right hemi
rh_highPGS = rh_highPGS.sum(axis=0)
diff_PGS_rh = rh_highPGS - rh_lowPGS
# same thing without zeros
lh_lowPGS_noNAN = lh_lowPGS_noNAN.sum(axis=0)
lh_highPGS_noNAN = lh_highPGS_noNAN.sum(axis=0)
diff_PGS_lh_noNAN = lh_highPGS_noNAN - lh_lowPGS_noNAN
rh_lowPGS_noNAN = rh_lowPGS_noNAN.sum(axis=0)  # right hemi
rh_highPGS_noNAN = rh_highPGS_noNAN.sum(axis=0)
diff_PGS_rh_noNAN = rh_highPGS_noNAN - rh_lowPGS_noNAN

# making random subsets

if 1 == 2:
    #holding_df = np.zeros((4, 10000))

    holding_df = np.zeros((163842, 100))

    for i in range(100):
        list_of_random_subs_1 = random.sample(range(550), 85)
    # need to make sure I am sampling without replacement
        new_list_to_pick_from = [a for a in range(550) if a not in list_of_random_subs_1]
        list_of_random_subs_2 = random.sample(new_list_to_pick_from, 85)
        df_rand1 = lh[list_of_random_subs_1]
        df_rand2 = lh[list_of_random_subs_2]

    # this could be in a fucked up subject order yet i don't think it matters cause we will take the mean
        df_rand1 = df_rand1.sum(axis=0)
        df_rand2 = df_rand2.sum(axis=0)

        df_delta = df_rand2 - df_rand1

        holding_df[:,i] = df_delta

        holding_df[:,i] = [df_delta.min(), df_delta.max(), df_delta.mean(), df_delta.std()]

    # I want the min, max, mean and SD of each permutation
    holding_df[0, :].mean() # min = -4.3
    holding_df[1, :].mean() # max = 4.3
    holding_df[2, :].mean() # mean = 0
    holding_df[3, :].mean() # sd = 1.17

# now you can make your color scale based on the min max from these and choose a representive random noise sample

# so I want to find out the right seed's to get min/max between 4-4.5, mean +/- .4 of 0 and std between 1-1.5

# this is a while loop, it will output the seed needed


if 42 == 231:
    logic = 0
    seedcount = 0
    while logic == 0:
        random.seed(seedcount)
        list_of_random_subs_1 = random.sample(range(550), 85)
        random.seed(seedcount)
        new_list_to_pick_from = [a for a in range(550) if a not in list_of_random_subs_1]
        list_of_random_subs_2 = random.sample(new_list_to_pick_from, 85)

        # this could be in a fucked up subject order yet I don't think it matters cause we will mean
        df_rand1 = lh[list_of_random_subs_1]
        df_rand2 = lh[list_of_random_subs_2]

        # getting the avg
        df_rand1 = df_rand1.sum(axis=0)
        df_rand2 = df_rand2.sum(axis=0)

        diff_random = df_rand2 - df_rand1

        logical_mean = -.2 < diff_random.mean() < .2
        logical_sd = 1 < diff_random.std() < 1.5

        if logical_mean == True & logical_sd == True:
            print(seedcount)  # its 5! lucky fuck
            logic = 1

        seedcount = seedcount + 1

# found random seed 5 to work
random.seed(5)
list_of_random_subs_1 = random.sample(range(550), 85)
random.seed(5)
new_list_to_pick_from = [a for a in range(550) if a not in list_of_random_subs_1]
list_of_random_subs_2 = random.sample(new_list_to_pick_from, 85)

# this could be in a fucked up subject order yet I don't think it matters cause we will mean
df_rand1 = lh[list_of_random_subs_1]
df_rand2 = lh[list_of_random_subs_2]
# same thing without zeros
df_rand1_noNAN = lh_noNAN[list_of_random_subs_1]
df_rand2_noNAN = lh_noNAN[list_of_random_subs_2]

# getting the avg
df_rand1 = df_rand1.sum(axis=0)
df_rand2 = df_rand2.sum(axis=0)
diff_random = df_rand2 - df_rand1
# same thing without zeros
df_rand1_noNAN = df_rand1_noNAN.sum(axis=0)
df_rand2_noNAN = df_rand2_noNAN.sum(axis=0)
diff_random_noNAN = df_rand2_noNAN - df_rand1_noNAN


# COORDINATE SECTION

# Import coordinates from sphere.reg
vertex_coords = fsio.read_geometry("/Users/nicjud/Projects/imaging/imagen/nick/analysis/py_vectorSubtraction/lh.sphere.reg")[0]

# Extract cartesian coordinates
x = vertex_coords[:, 0]
y = vertex_coords[:, 1]
z = vertex_coords[:, 2]

# Transform to spherical coordinates, radius is constant (100) for all points
theta = np.arccos(np.divide(z, 100))
phi = np.arctan2(y, x)

# PLOTTING SECTION

if 2 == 42:
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    # fig.suptitle('Vertex-wise global effects')
    # plot SES data above and below one SD
    ax1.tricontourf(theta, phi, diff_random, 15, vmin=-4.3, vmax=4.3)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    # now pgs above and below 1sd
    ax2.tricontourf(theta, phi, diff_PGS_lh, 15, vmin=-4.3, vmax=4.3)
    ax2.axes.get_xaxis().set_visible(False)
    ax2.axes.get_yaxis().set_visible(False)
    # now random data
    ax3.tricontourf(theta, phi, diff_SES_lh, 15, vmin=-4.3, vmax=4.3)
    ax3.axes.get_xaxis().set_visible(False)
    ax3.axes.get_yaxis().set_visible(False)
    ct = ax1.tricontourf(theta, phi, diff_random, 15)
    fig.colorbar(ct)
    fig.clim(-4.3, 4.3)
    plt.show()


# https://jdhao.github.io/2017/06/11/mpl_multiplot_one_colorbar/
from mpl_toolkits.axes_grid1 import AxesGrid
fig = plt.figure(figsize=(7, 4))

grid = AxesGrid(fig, 111,
                nrows_ncols=(1, 3),
                axes_pad=0.05,
                cbar_mode='single',
                cbar_location='left',
                cbar_pad=0.1
                )

colbar = ax.tricontourf(triang, diff_random, 15, vmin=-4.3, vmax=4.3)  # just for the color bar
colbar.set_cmap('plasma')
data_list = [diff_random, diff_PGS_lh, diff_SES_lh]
c = 0

for ax in grid:
    ax.set_axis_off()
    triang = tri.Triangulation(theta, phi)
    mask = np.all(np.where(lh_mask[triang.triangles], True, False), axis=1)
    triang.set_mask(mask)
    colplt = ax.tricontourf(triang, data_list[c], 15, vmin=-4.3, vmax=4.3)
    colplt.set_cmap('plasma')

    c = c + 1


cbar = ax.cax.colorbar(colbar)
#cbar = grid.cbar_axes[0].colorbar(colbar)
cbar.set_label_text("Surface Area")
#cbar.ax.set_yticklabels(['-4.3', '0', '4.3'])

#plt.savefig('/Users/nicjud/Projects/imaging/imagen/nick/analysis/figures/publication_graphs/flat_map_MAX.png', dpi=200)

# save in 1000 dpi

#cbar.ax.set_ybound(-4.3, 4.3)
#cbar.ax.set_ylim(-4.3, 4.3)




# density histogram
plt.show()
dens = sns.distplot(diff_random_noNAN,  kde=True, hist = False, label='random', color='red')
dens.figure.set_size_inches(10, 3)
dens.axes.get_yaxis().set_visible(False)
sns.distplot(diff_PGS_lh_noNAN,  kde=True, hist = False, label='PGS', color='green')
sns.distplot(diff_SES_lh_noNAN,  kde=True, hist = False, label='SES', color='blue')
sns.despine(top=True, left=True)
dens._remove_legend(True)
plt.savefig('/Users/nicjud/Projects/imaging/imagen/nick/analysis/figures/publication_graphs/dens_plot.png', dpi=200)

dens = sns.distplot(diff_random_noNAN,  kde=True, hist = False, label='random', color='red')
dens.figure.set_size_inches(10, 3)
dens.axes.get_yaxis().set_visible(False)
sns.distplot(diff_PGS_lh_noNAN,  kde=True, hist = False, label='PGS', color='green')
sns.distplot(diff_SES_lh_noNAN,  kde=True, hist = False, label='SES', color='blue')
sns.despine(top=True, left=True)
plt.savefig('/Users/nicjud/Projects/imaging/imagen/nick/analysis/figures/publication_graphs/dens_plot_legend.png', dpi=400)



########## playspace

# plotting on the brain

brain = Brain(subject_id='fsaverage', hemi='lh', surf='inflated')


#brain = Brain('fsaverage', 'inflated', hemi='lh',  views=['lat', 'med'])
brain.add_data(diff_SES,0, 8, center=0)

brain = Brain(subject_id='fsaverage', hemi='lh', surf='inflated', views=['lat', 'med'])
brain.add_data(diff_real,0, 8, center=0)


# now added a beta file with SES no global

SES_bl = mgh.load("/Users/nicjud/Projects/imaging/imagen/nick/analysis/py_vectorSubtraction/data/SES_TP1_Surf_lh_NOGLOBAL_beta.mgh")
SES_bl = SES_bl.get_data()

PGS_bl_nozeros = SES_bl[:,0,0,4]

SES_bl_nozeros = SES_bl[:,0,0,5]

sns.distplot(SES_bl[:,0,0,2],  kde=True, hist = False, label='Global', color='red')
sns.distplot(SES_bl[:,0,0,4],  kde=True, hist = False, label='PGS', color='green')
sns.distplot(SES_bl[:,0,0,5],  kde=True, hist = False, label='SES', color='blue')
plt.show()


# not needed since the PGS info is in the dimension of SES
#PGS_b = mgh.load("/Users/nicjud/Projects/imaging/imagen/nick/analysis/py_vectorSubtraction/data/PGS_TP1_Surf_lh_NOGLOBAL_beta.mgh")
#PGS_b = PGS_b.get_data()

PGS_bl = SES_bl[:,0,0,4]
SES_bl = SES_bl[:,0,0,5] # last is ses

brain = Brain(subject_id='fsaverage', hemi='lh', surf='inflated')
brain.add_data(PGS_bl, center=0)
brain = Brain(subject_id='fsaverage', hemi='lh', surf='inflated')
brain.add_data(SES_bl, center=0)

#from nilearn import plotting
#from nilearn import datasets
#fsaverage = datasets.fetch_surf_fsaverage5()


SES_br = mgh.load("/Users/nicjud/Projects/imaging/imagen/nick/analysis/py_vectorSubtraction/SES_TP1_Surf_rh_NOGLOBAL_beta.mgh")
SES_br = SES_br.get_data()



#plotting.plot_surf_stat_map(fsaverage.infl_right, diff, hemi='right',
#                            title='Surface right hemisphere', colorbar=True)


#brain = Brain("fsaverage", "left", "pial", views="frontal", background="dimgray")

"""
Because the morphometry files generated by
recon-all live in a predicatble location,
all you need to call the add_morphometry
method with is the name of the measure you want.
Here, we'll look at cortical curvatuve values,
and plot them for both hemispheres.
"""


"""
Each of the possible values is displayed in an
appropriate full-color map, but you can also
display in grayscale. Here we only plot the
left hemisphere.
"""
# brain.add_morphometry(diff, hemi='lh', grayscale=True)

"""
You can also use a custom colormap and tweak its range.
"""
#brain.add_morphometry("thickness",
#                      colormap="PuBuGn", min=1, max=4)