# Taking hemelb as base image
FROM hemelb
MAINTAINER Miguel O. Bernabeu (miguel.bernabeu@ed.ac.uk)

##
# Dependencies
##
RUN apt-get update && \
    apt-get install -y qhull-bin unzip wget libfreetype6-dev pkg-config python-tk && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN pip install matplotlib

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
# MCR installs its own linux runtime libraries and uses $LD_LIBRARY_PATH to give them priority over system ones.
# Python code running inside PolNet and trying to load the expat package finds the wrong version of libexpat.so.1 and fails at runtime.
##
RUN ln -fs /lib/x86_64-linux-gnu/libexpat.so.1 /opt/mcr/v85/bin/glnxa64/libexpat.so.1

##
# Download and install the standalone version of PolNet
##
WORKDIR /tmp
RUN wget https://www.dropbox.com/s/mn656wakrqkoiaz/PolNet_files.zip?dl=0 &&\
    mv PolNet_files.zip?dl=0 PolNet_files.zip && \
    unzip PolNet_files.zip && \
    cp PolNet_files/* /usr/local/bin/ && \
    chmod +x /usr/local/bin/PolNet /usr/local/bin/run_PolNet.sh && \
    rm -rf PolNet_files*

##
# Check out hemelb-hoff repo somewhere visible
##
WORKDIR /usr/local/bin
RUN git clone git@github.com:EPCCed/hemelb-hoff.git

##
# Place a PolNet launcher and a symlink to the data directory on the Desktop of the ubuntu user home space
##
RUN mkdir -p /home/ubuntu/Desktop/
COPY polnet.desktop /home/ubuntu/Desktop/
RUN ln -s /data /home/ubuntu/Desktop/

##
# Download and place the example data distributed with the protocol paper in a Desktop subdirectory of the ubuntu user home space
##
WORKDIR /tmp
RUN wget https://www.dropbox.com/s/nc8l6xig26jbw0a/990_Example2-skeleton.tif?dl=0 https://www.dropbox.com/s/n0pforgr7r7dc9p/990_Example2-flat.tif?dl=0 && \
    mv 990_Example2-skeleton.tif?dl=0 990_Example2-skeleton.tif && \
    mv 990_Example2-flat.tif?dl=0 990_Example2-flat.tif && \
    mkdir /home/ubuntu/Desktop/paper_example_data/ && \
    mv 990_Example2* /home/ubuntu/Desktop/paper_example_data/ && \
    echo "Any intermediate file saved to this directory will disappear every time that the container is stopped and restarted. Copy the example data to /data (which maps to a directory in the host machine) in order for the intermediate files to persist." > /home/ubuntu/Desktop/paper_example_data/README.txt

# The MATLAB runtime overwrites $LD_LIBRARY_PATH ignoring its previous content. Copy VTK/VMTK libraries to a visible location
RUN cp $VMTKHOME/lib/lib* /usr/lib/x86_64-linux-gnu/

##
# Run PolNet at login
##
#RUN mkdir /etc/xdg/autostart/
#COPY polnet.desktop /etc/xdg/autostart/
