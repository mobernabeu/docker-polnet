# Taking hemelb as base image
FROM hemelb
MAINTAINER Miguel O. Bernabeu (miguel.bernabeu@ed.ac.uk)

##
# Dependencies
##
RUN apt-get update && \
    apt-get install -y qhull-bin unzip wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

##
# Download and install MATLAB's MCR
##
WORKDIR /opt
RUN wget http://www.mathworks.com/supportfiles/downloads/R2015a/deployment_files/R2015a/installers/glnxa64/MCR_R2015a_glnxa64_installer.zip && \
    unzip MCR_R2015a_glnxa64_installer.zip && \
    mkdir /opt/mcr && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    rm MCR_R2015a_glnxa64_installer.zip

##
# Download and install the standalone version of PolNet
##
WORKDIR /tmp
RUN wget https://www.dropbox.com/s/ru2l14aledjho6r/PolNet_files.zip?dl=0 && \
    mv PolNet_files.zip?dl=0 PolNet_files.zip && \
    unzip PolNet_files.zip && \
    cp PolNet_files/* /usr/local/bin/ && \
    chmod +x /usr/local/bin/PolNet /usr/local/bin/run_PolNet.sh && \
    rm -rf PolNet_files*

##
# Place a PolNet launcher and a symlink to the data directory on the Desktop of any new user
##
RUN mkdir /etc/skel/Desktop/
COPY polnet.desktop /etc/skel/Desktop/
RUN ln -s /data /etc/skel/Desktop/

##
# Run PolNet at login
##
#RUN mkdir /etc/xdg/autostart/
#COPY polnet.desktop /etc/xdg/autostart/
