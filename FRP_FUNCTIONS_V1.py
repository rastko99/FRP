#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def linear(x,a,b):
    """Models a linear function"""
    return (a*x + b)

def RFP(uncal_file,output_dir,output_name,allframes = False ,stages = 3):
    """
    Runs Full Pipeline specified for MIRI
    
    Input:
    -uncal_file: Path to the uncalibrated datafile.
    -output_dir: Path to the desired output directory.
    -output_name: name for the output data file.
    -allframes: Boolean input parameter,
                if True all of the frames available are used by skipping the first- and lastframe steps,
                if False the first and last frames are omitted as per STScI recommendation.
    -stages: Integer value specifying up to which stage the pipeline should be ran.
             E.g. if stages = 2 only runs the pipeline up to and including stage2.
             
    Output:
    -Returns the output datafile with name equal to output_name+suffix in output_dir corresponding to what stages have been run.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
            
    detector1 = calwebb_detector1.Detector1Pipeline()

    detector1.output_dir = output_dir
    detector1.save_results = True
    detector1.output_file = output_name
    detector1.jump.use_ellipses = True  
    if allframes:
        detector1.firstframe.skip = True
        detector1.lastframe.skip = True

    detector1_output = detector1.run(uncal_file)
    
    if stages >= 2:
        
        image2 = calwebb_image2.Image2Pipeline()

        image2.output_dir = output_dir
        image2.save_results = True
        image2.output_file = output_name

        image2_output = image2.run(detector1_output)
    
    if stages >= 3:
        
        image3 = calwebb_image3.Image3Pipeline()

        image3.output_dir = output_dir
        image3.save_results = True
        image3.output_file = output_name

        image3.run(image2_output)
    return

def BackgroundNoise(filepath, windowsize = [10,24,50] ,plot = False):
    """
    Calculates Noise (stdev) and Background (average) values for 6 regions (splices) of the MIRI detector.
    
    Input: 
    -filepath: A filepath leading to the input .fits file.
    -windowsize: List of even windowsizes over which to comptue background and noise, values specified should not be higher than 50.
    -plot: Boolean variable, if True plots the windows over which the background and noise have been calculated, along with the calculated values
    
    Output:
    -Returns tuple of 2D-arrays containing background and noise values per windowsize.
     Individual arrays are of shape (3,7) each row being the values for each windowsize,
     each column being the values per region,
     the 7th column being the average of all those regions.
    """
    
    #Defines region centers for regions 1-6.
    
    center_x = np.array([400,873,560,790,605,815])
    center_y = np.array([900,905,485,600,120,200])
        
    #Data read-in and manipulation
    data =  fits.getdata(filepath,"SCI")
    flat_data = data.flatten()
    background = np.zeros(0)
    noise = np.zeros(0)
    for window in windowsize:
        
        x_start, y_start = (center_x - window/2).astype(int), (center_y - window/2).astype(int)
        x_end, y_end = (center_x + window/2).astype(int), (center_y + window/2).astype(int)
        
        region_1 = data[y_start[0]:y_end[0],x_start[0]:x_end[0]]
        region_2 = data[y_start[1]:y_end[1],x_start[1]:x_end[1]]
        region_3 = data[y_start[2]:y_end[2],x_start[2]:x_end[2]]
        region_4 = data[y_start[3]:y_end[3],x_start[3]:x_end[3]]
        region_5 = data[y_start[4]:y_end[4],x_start[4]:x_end[4]]
        region_6 = data[y_start[5]:y_end[5],x_start[5]:x_end[5]]
        
        background_1, noise_1 = np.average(region_1.flatten()), np.std(region_1.flatten())
        background_2, noise_2 = np.average(region_2.flatten()), np.std(region_2.flatten())
        background_3, noise_3 = np.average(region_3.flatten()), np.std(region_3.flatten())
        background_4, noise_4 = np.average(region_4.flatten()), np.std(region_4.flatten())
        background_5, noise_5 = np.average(region_5.flatten()), np.std(region_5.flatten())
        background_6, noise_6 = np.average(region_6.flatten()), np.std(region_6.flatten())
        
        background_regions = np.array([background_1,background_2,background_3,background_4,background_5,background_6])
        background = np.append(background,background_regions)
        
        noise_regions = np.array([noise_1,noise_2,noise_3,noise_4,noise_5,noise_6])
        noise = np.append(noise,noise_regions)
        
        background, noise = np.append(background, np.average(background_regions)), np.append(noise, np.average(noise_regions))

        if plot:
            fig,[[ax1,ax2],[ax3,ax4],[ax5,ax6]] = plt.subplots(3,2,sharex = True,sharey= True,figsize = (10,10), layout = "tight")
            
            fig.suptitle(f"{window} pixel window",fontsize = 18)
            ax1.imshow(region_1)
            ax1.set_title(f"Region 1 \n Background : {background_1:0.3f} Noise : {noise_1:0.3f} ")
            ax2.imshow(region_2)
            ax2.set_title(f"Region 2 \n Background : {background_2:0.3f} Noise : {noise_2:0.3f} ")
            ax3.imshow(region_3)
            ax3.set_title(f"Region 3 \n Background : {background_3:0.3f} Noise : {noise_3:0.3f} ")
            ax4.imshow(region_4)
            ax4.set_title(f"Region 4 \n Background : {background_4:0.3f} Noise : {noise_4:0.3f} ")
            ax5.imshow(region_5)
            ax5.set_title(f"Region 5 \n Background : {background_5:0.3f} Noise : {noise_5:0.3f} ")
            ax6.imshow(region_6)
            ax6.set_title(f"Region 6 \n Background : {background_6:0.3f} Noise : {noise_6:0.3f} ")
    
    background, noise = background.reshape(3,7), noise.reshape(3,7)
    return background,noise

def dithercombine(asn_file,method,output_dir,output_name):
    """
    Combines dithers according to the information in the provided association file using the appropriate steps in the stage3 pipeline.
    
    Input:
    -asn_file: Path the to association file containing relevant information on the dithers.
    -method: Parameter specifying the skymatching method to be used when combining the dithers and subtracting the sky.
    -output_dir: Path to the desired output directory.
    -output_name: Name for the output data file.
    
    Output:
    -Returns the dithercombined result of the stage3 pipeline in the output_dir named output_name+suffix.
    """
    image3 = calwebb_image3.Image3Pipeline()

    image3.output_dir = output_dir
    image3.save_results = True
    image3.output_file = output_name
    image3.skymatch.subtract = True
    image3.skymatch.skymethod = method
    image3.skymatch.skystat ="mean"
    #image3.skymatch.match_down = False #default is True
    #image3.tweakreg.skip = True #this had to be added for medsubbed cal files to be asn combined
    #image3.source_catalog.skip = True #this had to be added for medsubbed cal files to be asn combined
    #image3.skymatch.skip = True
    #image3.outlier_detection.skip = True
    
    
    image3.run(asn_file)
    return

def Noise_2(filepath, kernelsize = 5):
    """
    Second iteration of the Noise Estimator. 
    Estimates the noise using by looping a kernel over three preselected regions.
    This effectively increases the number of subregions significantly
    
    Input:
    -filepath: Path to the file for which the noise should be calculated
    -kernelsize: Integer value of which kernelsize will be used, only accepts kernelsizes 3,4,5 and multiplicatives of those values.
    
    Output:
    -Returns a total array containing three region arrays.
     Each region array within the total array contains the separate noise values calculated on each kernel.
     Making conclusive statements on the noise values is more easily done by taking the average of each of the region arrays.
    
    """
    data = fits.getdata(filepath, "SCI")
    region_1 = data[130:250,770:950]
    region_2 = data[688:808,810:930]
    region_3 = data[60:120,880:1000]
    kernel = kernelsize
    noise_arr_1,noise_arr_2,noise_arr_3 = np.zeros(0),np.zeros(0),np.zeros(0)
    
    for i in range(int(len(region_1[:,0])/kernel)):
        for j in range(int(len(region_1[0,:])/kernel)):
            subregion = region_1[i*kernel:(i+1)*kernel,j*kernel:(j+1)*kernel]
            noise = np.std(subregion)
            noise_arr_1 = np.append(noise_arr_1,noise)
    for i in range(int(len(region_2[:,0])/kernel)):
        for j in range(int(len(region_2[0,:])/kernel)):
            subregion = region_2[i*kernel:(i+1)*kernel,j*kernel:(j+1)*kernel]
            noise = np.std(subregion)
            noise_arr_2 = np.append(noise_arr_2,noise)
    for i in range(int(len(region_3[:,0])/kernel)):
        for j in range(int(len(region_3[0,:])/kernel)):
            subregion = region_3[i*kernel:(i+1)*kernel,j*kernel:(j+1)*kernel]
            noise = np.std(subregion)
            noise_arr_3 = np.append(noise_arr_3,noise)
    noise_arr = ([noise_arr_1,noise_arr_2,noise_arr_3])
    
    return noise_arr

def bkg_sub(filepath,kernel = 5):
    """
    Performs manual background subtraction based on lower right region
    
    Input:
    -filepath: Path to the datafile for which the background should be subtracted.
    -kernel: Integer value used for the kernels to estimate the background and perform subtraction.
    
    Output:
    -Returns the background subtracted datafile in the working directory
    """
    
    if "medsub" in filepath:
        print("Check file(name) probably is already bkg_subbed")
        return
    data = fits.getdata(filepath,"SCI")
    region_3 = data[60:120,880:1000]
    bkg_array = np.zeros(0)
    for i in range(int(len(region_3[:,0])/kernel)):
        for j in range(int(len(region_3[0,:])/kernel)):
            subregion = region_3[i*kernel:(i+1)*kernel,j*kernel:(j+1)*kernel]
            bkg = np.median(subregion)
            bkg_array = np.append(bkg_array,bkg)
    bkg_subbed_data = data - np.median(bkg_array)
    shutil.copyfile(filepath,filepath.replace("cal","medsub_cal"))
    newfile = fits.open(filepath.replace("cal","medsub_cal"),mode = "update")
    newfile["SCI"].data = bkg_subbed_data
    newfile.close()

    return print(np.std(bkg_array))

def all_CMD_plotter(red,green,blue,savefig = ""):
    """
    Plots all permutations of the CMD given three filters.
    
    Input:
    -red: Path to the largest wavelength filter data (F1500W).
    -green: Path to the middle wavelength filter data (F1000W).
    -blue: Path to the smallest wavelength filter data (F770W).
    -savefig: String containing the name of the figure to be saved. If left empty the figure is not saved.
    
    Output:
    -Plots all permutations of the CMD given the three filters, if savefig is specified the figure is also saved to the working directory.
    """
    blue_fits,red_fits,green_fits = fits.open(blue),fits.open(red),fits.open(green)
    blue_data,red_data,green_data = blue_fits[1].data,red_fits[1].data,green_fits[1].data
    blue_filter,red_filter,green_filter = str((blue_data.columns[-2]).name),str((red_data.columns[-2]).name),str((green_data.columns[-2]).name)
    blue_mags,red_mags,green_mags = blue_data[blue_filter],red_data[red_filter],green_data[green_filter]
    
    fig,axs = plt.subplots(2,3,sharex= "col",figsize = (16,8))
    fig.suptitle("Color Magnitude Diagrams")
    #fig.tight_layout()
    
    axs[0,0].scatter(blue_mags-red_mags,red_mags,alpha = 0.5,edgecolor= "black")
    axs[0,0].grid(alpha=0.5)
    axs[0,0].set_ylabel(f"{red_filter}")    
    axs[0,0].set_ylim(axs[0,0].get_ylim()[::-1])
    axs[0,0].set_title(f"{blue_filter}-{red_filter}")
    
    axs[1,0].scatter(blue_mags-red_mags,blue_mags,alpha = 0.5,edgecolor= "black")
    axs[1,0].grid(alpha=0.5)
    axs[1,0].set_xlabel(f"{blue_filter}-{red_filter}")
    axs[1,0].set_ylabel(f"{blue_filter}")
    axs[1,0].set_ylim(axs[1,0].get_ylim()[::-1])
    
    axs[0,1].scatter(green_mags-red_mags,red_mags,alpha = 0.5,edgecolor= "black")
    axs[0,1].grid(alpha=0.5)
    axs[0,1].set_ylabel(f"{red_filter}")
    axs[0,1].set_ylim(axs[0,1].get_ylim()[::-1])
    axs[0,1].set_title(f"{green_filter}-{red_filter}")
    
    axs[1,1].scatter(green_mags-red_mags,green_mags,alpha = 0.5,edgecolor= "black")
    axs[1,1].grid(alpha=0.5)
    axs[1,1].set_xlabel(f"{green_filter}-{red_filter}")
    axs[1,1].set_ylabel(f"{green_filter}") 
    axs[1,1].set_ylim(axs[1,1].get_ylim()[::-1])
    
    axs[0,2].scatter(blue_mags-green_mags,green_mags,alpha = 0.5,edgecolor= "black")
    axs[0,2].grid(alpha=0.5)
    axs[0,2].set_ylabel(f"{green_filter}")    
    axs[0,2].set_ylim(axs[0,2].get_ylim()[::-1])
    axs[0,2].set_title(f"{blue_filter}-{green_filter}")
    
    axs[1,2].scatter(blue_mags-green_mags,blue_mags,alpha = 0.5,edgecolor= "black")
    axs[1,2].grid(alpha=0.5)
    axs[1,2].set_xlabel(f"{blue_filter}-{green_filter}")
    axs[1,2].set_ylabel(f"{blue_filter}") 
    axs[1,2].set_ylim(axs[1,2].get_ylim()[::-1])
    
    fig.tight_layout()
    
    if savefig != "":
        fig.savefig(savefig)

    blue_fits.close
    red_fits.close
    green_fits.close
    return

def CMD(blue,red):
    """
    Plots a CMD given two filters.
    
    Input:
    -blue: Filepath to the blue (smaller wavelength) FITS data.
    -red: Filepath to the red (larger wavelength) FITS data.
     
    Output:
    -Matplotlib Colour Magnitude Diagram.
    """
    blue_fits,red_fits = fits.open(blue),fits.open(red)
    blue_data,red_data = blue_fits[1].data,red_fits[1].data
    blue_filter,red_filter = str((blue_data.columns[-2]).name),str((red_data.columns[-2]).name)
    blue_mags,red_mags = blue_data[blue_filter],red_data[red_filter]
    
    fig,ax = plt.subplots()
    ax.scatter(blue_mags-red_mags,red_mags,alpha = 0.5,edgecolor= "black")
    ax.grid(alpha=0.5)
    ax.set_title("Color Magnitude Diagram")
    ax.set_xlabel(f"{blue_filter}-{red_filter}")
    ax.set_ylabel(f"{red_filter}")
    ax.set_ylim(ax.get_ylim()[::-1])
    
    blue_fits.close()
    red_fits.close()
    return

def sharpround(fitsfile,sbins = 50, rbins = 50):
    """
    Evaluatory tool that plots the sharpness and roundness of the detected sources.
    
    Input:
    -fitsfile: Filepath to FITS table containing detections made by starbug2.
    -sbins: Number of bins to use in plotting the roundness histogram.
    -rbins: Number of bins to use in plotting the sharpness histogram.
    
    Output:
    -Histogram of sharpness
    -Histogram of roundness
    """
    fitsdata = fits.open(fitsfile)
    data = fitsdata[1].data
    sharpness = data["sharpness"]
    r1,r2 = data["roundness1"],data["roundness2"]
    
    fig, (ax1,ax2) = plt.subplots(1,2,figsize = (16,6))
    
    ax1.hist(sharpness,bins =sbins)
    ax1.set_title("Sharpness")
    ax1.grid(alpha=0.5)
    
    ax2.hist(r1,color = "black",bins= rbins,label = "Roundness1")
    ax2.hist(r2,color = "red" ,bins = rbins, alpha = 0.5,label="Roundness2")
    ax2.set_title("Roundness")
    ax2.grid(alpha=0.5)
    ax2.legend()
    
    return

def regionmaker(F770W,F1000W,F1500W,outputname, threshold =0.0001):
    """
    Makes a region file that consists of sources that appear in at least 2 filters.
    
    Inputs:
    -F770W: Filepath to region list (.reg) corresponding to detected sources in F770W.
    -F1000W: Filepath to region list (.reg) corresponding to detected sources in F1000W.
    -F1500W: Filepath to region list (.reg) corresponding to detected sources in F1500W.
    -outputname: Outputname of the final region list.
    -threshold: threshold value in arcseconds within which the sources are counted as the same.
    
    Outputs:
    -Verbose output of how many of the total sources are omitted.
    -One final matched .reg region file.
    """
    reg77,reg10,reg15 = np.loadtxt(F770W,dtype = str,skiprows =1),np.loadtxt(F1000W,dtype = str,skiprows =1),np.loadtxt(F1500W,dtype = str,skiprows =1)
        
    ra_77,ra_10,ra_15 = reg77[:,1].astype(float),reg10[:,1].astype(float),reg15[:,1].astype(float)
    dec_77,dec_10,dec_15 = reg77[:,2].astype(float),reg10[:,2].astype(float),reg15[:,2].astype(float)
    
    size_77,size_10,size_15 = np.zeros(len(dec_77)),np.zeros(len(dec_10)),np.zeros(len(dec_15))
    
    for i in range(len(size_77)):
        size_77[i] = float(reg77[i,3].replace("i",""))
    for i in range(len(size_10)):
        size_10[i] = float(reg10[i,3].replace("i",""))
    for i in range(len(size_15)):
        size_15[i] = float(reg15[i,3].replace("i",""))
        
    ra_,dec_,size_ = np.zeros(0),np.zeros(0),np.zeros(0)
    size_ = size_.astype(str)
    for i in range(len(ra_77)):
        
        if np.logical_or(np.logical_and(np.any(np.abs(ra_10-ra_77[i]) <= (threshold*2)), np.any(np.abs(dec_10-dec_77[i]) <= threshold)),np.logical_and(np.any(np.abs(ra_15-ra_77[i]) <= (threshold*2)), np.any(np.abs(dec_15-dec_77[i]) <= threshold)) ):
            ra_ = np.append(ra_,ra_77[i])
            dec_ = np.append(dec_,dec_77[i])
            size_ = np.append(size_,str(size_77[i])+"i")
    
    firstcol = np.full(len(ra_), "fk5;circle")
    reg_final = np.column_stack((firstcol,ra_,dec_,size_))
    
    np.savetxt(fname = outputname+".reg", X = reg_final,fmt = "%0s", delimiter=' ', header = "global color=red width=2",comments = "")
    print(f"Out of a total of {len(ra_77)} sources in 7.7 micron, {(len(ra_77)-len(ra_))} / {((len(ra_77)-len(ra_)))*100/len(ra_77):.3f} % were omitted resulting in a total of {len(ra_)} sources left in the final region.")
    return 

def all_CMD_plotter_2(data,savefig = ""):
    """
    Plots all permutations of the CMD given three filters.
    
    Input:
    -data: FITS table containing aperture photometry of all three filters, as is obtained from starbug2 matching the seperate files.
    -savefig: String containing the name of the figure to be saved. If left empty the figure is not saved.
    
    Output:
    -Plots all permutations of the CMD given the three filters, if savefig is specified the figure is also saved to the working directory.
    """
    fitsdata = fits.open(data)
    fitsdata = fitsdata[1].data
    blue_filter,red_filter,green_filter = str((fitsdata.columns[-6]).name),str((fitsdata.columns[-2]).name),str((fitsdata.columns[-4]).name)
    blue_mags,red_mags,green_mags = fitsdata[blue_filter],fitsdata[red_filter],fitsdata[green_filter]
    
    fig,axs = plt.subplots(2,3,sharex= "col",figsize = (16,8))
    fig.suptitle("Color Magnitude Diagrams")
    #fig.tight_layout()
    
    axs[0,0].scatter(blue_mags-red_mags,red_mags,alpha = 0.5,edgecolor= "black")
    axs[0,0].grid(alpha=0.5)
    axs[0,0].set_ylabel(f"{red_filter}")    
    axs[0,0].set_ylim(axs[0,0].get_ylim()[::-1])
    axs[0,0].set_title(f"{blue_filter}-{red_filter}")
    
    axs[1,0].scatter(blue_mags-red_mags,blue_mags,alpha = 0.5,edgecolor= "black")
    axs[1,0].grid(alpha=0.5)
    axs[1,0].set_xlabel(f"{blue_filter}-{red_filter}")
    axs[1,0].set_ylabel(f"{blue_filter}")
    axs[1,0].set_ylim(axs[1,0].get_ylim()[::-1])
    
    axs[0,1].scatter(green_mags-red_mags,red_mags,alpha = 0.5,edgecolor= "black")
    axs[0,1].grid(alpha=0.5)
    axs[0,1].set_ylabel(f"{red_filter}")
    axs[0,1].set_ylim(axs[0,1].get_ylim()[::-1])
    axs[0,1].set_title(f"{green_filter}-{red_filter}")
    
    axs[1,1].scatter(green_mags-red_mags,green_mags,alpha = 0.5,edgecolor= "black")
    axs[1,1].grid(alpha=0.5)
    axs[1,1].set_xlabel(f"{green_filter}-{red_filter}")
    axs[1,1].set_ylabel(f"{green_filter}") 
    axs[1,1].set_ylim(axs[1,1].get_ylim()[::-1])
    
    axs[0,2].scatter(blue_mags-green_mags,green_mags,alpha = 0.5,edgecolor= "black")
    axs[0,2].grid(alpha=0.5)
    axs[0,2].set_ylabel(f"{green_filter}")    
    axs[0,2].set_ylim(axs[0,2].get_ylim()[::-1])
    axs[0,2].set_title(f"{blue_filter}-{green_filter}")
    
    axs[1,2].scatter(blue_mags-green_mags,blue_mags,alpha = 0.5,edgecolor= "black")
    axs[1,2].grid(alpha=0.5)
    axs[1,2].set_xlabel(f"{blue_filter}-{green_filter}")
    axs[1,2].set_ylabel(f"{blue_filter}") 
    axs[1,2].set_ylim(axs[1,2].get_ylim()[::-1])
    
    fig.tight_layout()
    
    if savefig != "":
        fig.savefig(savefig)
        
    return

