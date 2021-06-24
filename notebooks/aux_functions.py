import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np

# import seaborn as sns; sns.set(color_codes=True)
# sns.set_style("darkgrid")

def plot_color_bin(z,mean,sigma,dz=0,scatter_mean=False,label='RS Model',lcolor='r',axs=None):

    b = np.arange(0.08,0.92,0.03)
    indices_list = list(chunks2(z, b))

    y_bin = np.array([ np.nanmedian(mean[idx]) for idx in indices_list])
    # std_bin = np.array([(np.nanmedian(sigma[idx])**2+np.nanstd(mean[idx])**2)**(1/2) for idx in indices_list])
    std_bin = np.array([np.nanmedian(sigma[idx]) for idx in indices_list])
    if scatter_mean:
        std_bin = np.array([ np.nanstd(mean[idx]) for idx in indices_list])

    x_bin = np.array([ np.median(z[idx]) for idx in indices_list])

    if axs is None:
        plt.errorbar(x_bin,y_bin,yerr=std_bin,color=lcolor,label=label)
    else:
        axs.errorbar(dz+x_bin,y_bin,yerr=std_bin,color=lcolor,label=label)


def validating_color_model_grid(cat,cat2,color_list,lcolor=None,sigma=False,fraction=False):
	if lcolor is None: lcolor = color_list

	plt.clf()
	fig = plt.figure(figsize=(16,10))
	fig.subplots_adjust(hspace=0.03, wspace=0.35)
	# fig.tight_layout()

	for i,li in enumerate(color_list):
		axs = fig.add_subplot(2, 3, i+1)
		
		# #### Color plots
		# color = gal[lcolor[i]]
		# pmem = gal['Pmem']

		zrs = cat['redshift']
		mur = cat['rs_param_%s'%(li)][:,0]
		sigr = cat['rs_param_%s'%(li)][:,1]
		alpr = cat['rs_param_%s'%(li)][:,2]

		mub = cat['bc_param_%s'%(li)][:,0]
		sigb = cat['bc_param_%s'%(li)][:,1]
		alpb = cat['bc_param_%s'%(li)][:,2]

		zcls2 = cat2['redshift']
		mur2 = cat2['rs_param_%s'%(li)][:,0]
		sigr2 = cat2['rs_param_%s'%(li)][:,1]
		alpr2 = cat2['rs_param_%s'%(li)][:,2]

		mub2 = cat2['bc_param_%s'%(li)][:,0]
		sigb2 = cat2['bc_param_%s'%(li)][:,1]
		alpb2 = cat2['bc_param_%s'%(li)][:,2]

		# scale = 10*pmem**(1/2)
		if (not sigma)&(not fraction):
			# axs.scatter(zcls,color,s=scale,color='k',alpha=0.08)
			axs.scatter(zrs,mur,color='#F1948A',alpha=0.4,s=20,label='RS: individual systems')
			axs.scatter(zrs,mub,color='#85C1E9',alpha=0.4,s=20,label='BC: individual systems')
			
			idx = np.argsort(zcls2)
			plt.plot(zcls2[idx],mur2[idx],color='r',linestyle='--',linewidth=3)
			plt.plot(zcls2[idx],mub2[idx],color='b',linestyle='--',linewidth=3)

			# plot_color_bin(zrs,mur,sigr,scatter_mean=True,lcolor='#F1948A')
			# plot_color_bin(zrs,mub,sigb,scatter_mean=True,label='blue cloud',lcolor='#85C1E9')
			# plt.xlim(0.08,0.92)
			axs.set_ylim(0.,np.max(mur)+1.2*np.max(sigr))
			#axs.set_xlim(0.05,0.95)
			axs.set_ylabel(r'$%s$'%(lcolor[i]),fontsize=16)

		if sigma:
			axs.scatter(zrs,sigr,color='#F1948A',alpha=0.4,s=20,label='RS: individual systems')
			axs.scatter(zrs,sigb,color='#85C1E9',alpha=0.4,s=20,label='BC: individual systems')

			idx = np.argsort(zcls2)
			plt.plot(zcls2[idx],sigr2[idx],color='r',linestyle='--',linewidth=3)
			plt.plot(zcls2[idx],sigb2[idx],color='b',linestyle='--',linewidth=3)

			# plot_color_bin(zrs,sigr,sigr,scatter_mean=True,lcolor='#F1948A')
			# plot_color_bin(zrs,sigb,sigb,scatter_mean=True,label='blue cloud',lcolor='#85C1E9')
			# plt.xlim(0.08,0.92)
			# axs.set_ylim(0.,np.max(sigr)+1.2*np.max(sigr))
			# axs.set_yscale('log')
			#axs.set_xlim(0.05,0.95)
			axs.set_ylabel(r'$\sigma_{%s}$'%(lcolor[i]),fontsize=16)


		if fraction:
			axs.scatter(zrs,alpr,color='#F1948A',alpha=0.4,s=20,label='RS: individual systems')
			axs.scatter(zrs,alpb,color='#85C1E9',alpha=0.4,s=20,label='BC: individual systems')

			idx = np.argsort(zcls2)
			plt.plot(zcls2[idx],alpr2[idx],color='r',linestyle='--',linewidth=3)
			plt.plot(zcls2[idx],alpb2[idx],color='b',linestyle='--',linewidth=3)

			# plot_color_bin(zrs,alpr,sigr,scatter_mean=True,lcolor='#F1948A')
			# plot_color_bin(zrs,alpb,sigb,scatter_mean=True,label='blue cloud',lcolor='#85C1E9')
			# plt.xlim(0.08,0.92)
			# axs.set_ylim(0.,np.max(sigr)+1.2*np.max(sigr))
			#axs.set_xlim(0.05,0.95)
			axs.set_ylabel(r'$w {%s}$'%(lcolor[i]),fontsize=16)

		if i>=3:
			axs.set_xlabel('redshift',fontsize=16)

def color_redshift_red_blue_bad(gal,color_label,lcol="color",delta=False):
	zcls = gal["redshift"]
	color= gal[lcol]
	pmem=gal["Pmem"]
	pz = gal["Pz"]

	w, = np.where((gal["red_score"]<=2)&((gal["rs_bad"])<=1))
	w2, = np.where((gal["red_score"]>2)&((gal["rs_bad"])<=1))
	w3, = np.where((gal["rs_bad"]>1))

	plt.clf()
	fig = plt.figure(figsize=(16,8))
	fig.subplots_adjust(hspace=0.03, wspace=0.35)
	# fig.tight_layout()

	for i,li in enumerate(color_label):
		c_max = np.mean(color[w2,i])+5*np.std(color[w2,i])
		c_min = np.mean(color[w,i])-3*np.std(color[w,i])

		axs = fig.add_subplot(2, 3, i+1)

		axs.scatter(zcls[w],color[w,i],color="b",alpha=0.1,s=10)
		axs.scatter(zcls[w2],color[w2,i],color="r",alpha=0.05,s=10)
		
		if not delta:
			axs.scatter(zcls[w3],color[w3,i],color="k",alpha=0.05,s=10)
			cmax = np.mean(color[w2,i])+3/5*np.std(color[w2,i])

		axs.set_ylim(c_min,c_max)
		
		label = r'$%s$'%(color_label[i])
		if delta: label = r'$\Delta %s$'%(color_label[i])
		axs.set_ylabel(label,fontsize=16)
		axs.set_xlim(0.05,0.95)

		if i>=3:
			axs.set_xlabel('redshift',fontsize=16)

def color_hist_red_blue(gal,color_label,lcol="color",delta=True):
	zcls = gal["redshift"]
	color= gal[lcol]
	pmem=gal["Pmem"]
	pz = gal["Pz"]
	ntotal=np.sum(pmem)

	w, = np.where((gal["red_score"]<=2)&((gal["rs_bad"])<=1))
	w2, = np.where((gal["red_score"]>2)&((gal["rs_bad"])<=1))
	w3, = np.where((gal["rs_bad"]>1))

	plt.clf()
	fig, axis = plt.subplots(2, 3, sharey=True,sharex=False, figsize=(14,10))
	fig.subplots_adjust(hspace=0.45, wspace=0.05)
	# fig.tight_layout()

	j=0
	for i,li in enumerate(color_label):
		i2=i
		if i>2:i2=i-3; j=1

		axs = axis[j][i2]

		c_max = np.mean(color[w2,i])+5*np.std(color[w2,i])
		c_min = np.mean(color[w,i])-3*np.std(color[w,i])

		nbins = np.linspace(c_min,c_max,25)

		axs.hist(color[w2,i],bins=nbins,color="r",alpha=1.,weights=pmem[w2]/ntotal)
		axs.hist(color[w,i],bins=nbins,color="b",alpha=0.7,weights=pmem[w]/ntotal)
		
		# if not delta:
		# axs.hist(color[w3,i],bins=nbins,color="k",alpha=0.9,weights=pmem[w3]/ntotal)
		
		label = r'$%s$'%(color_label[i])
		if delta: label = r'$\Delta %s$'%(color_label[i])
		axs.set_xlabel(label,fontsize=16)
		# axs.set_ylim(-0.05,ymax)
		axs.set_xlim(c_min,c_max)

		if (i==0)|(i==3):
			axs.set_ylabel('fraction',fontsize=16)
