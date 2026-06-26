# Use the pre-built tool stash image as the base
FROM parkinsonlab/mpro_stash:v1.0.8

# Download and install the DeepProZyme wrapper
RUN wget https://compsysbio.org/metapro_hosting/deepec_wrapper.py -O /usr/local/bin/run_deepectransformer.py \
&& chmod +x /usr/local/bin/run_deepectransformer.py

# Set the working directory for the pipeline code
WORKDIR /pipeline

ARG CACHE-BUSTER=2
RUN echo "busting cache here"

# Clone the pipeline code from the 'test' branch into the current directory
# This single RUN command is more efficient than dozens of ADD commands
RUN git clone --branch v4.0.2 https://github.com/ParkinsonLab/MetaPro.git .

# Set permissions for the pipeline scripts
RUN chmod -R 555 /pipeline

# Set the final working directory and default command
WORKDIR /
CMD ["bash"]