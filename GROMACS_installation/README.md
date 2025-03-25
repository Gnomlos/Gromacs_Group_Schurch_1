# GROMACS intallation
In this part of the repo, is the tutorial of how to intall GROMACS in linux. Be aware that the linux distribution used is Fedora. This implies that some of the command needs to be adjusted to the version of the linux used.
# Setting up the environement
The first think to do is to set-up the C and C++ environement of your system. If this is alreay done you can skip to the next point. If not then, you will need to install Cmake, gcc and g++

# 1. Download GROMACS and prepapring the files for the installation
The first thing you need to do is to download GROMACS. Go to https://manual.gromacs.org/2024.4/download.html and choose the wanted version. In our case we iinstall the 2024.4v. Once downloaded the you need to extract the files. This can be done in 2 ways. The first one is using you graphical interface to extract it to the desired location. The second one is through the console, using the command "tar xfz" followed by your file name.

Then you will navigate inside your file and create a directory "build" which you will get into. 

# 2. Intallation of GROOMACS
Now we will install GROMACS on our computer:
1. use the "cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON" with the path of the top directory of GROMACS in place of the 2 dots.
2. type "make"
3. type "make check"
4. type "sudo make install"
5. type "source /usr/local/gromacs/bin/GMXRC" <-- This command allows us to use GROMACS.

