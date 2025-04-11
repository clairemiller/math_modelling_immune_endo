FROM gcc:14.2.0

# Make python a symbolic link to python3
RUN ln -s /usr/bin/python3 /usr/bin/python

# Clone the github repo into /project/auto-07p
RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*
WORKDIR /project
RUN git clone https://github.com/auto-07p/auto-07p/

# Install auto-07 and add to the python path
WORKDIR /project/auto-07p
RUN ./configure --prefix=/project --disable-plaut04-qt
RUN make && make clean
ENV PYTHONPATH=$PYTHONPATH:/project/auto-07p/python

# Install the appropriate packages
RUN apt-get update && apt-get install -y python3-numpy python3-pandas

# Copy the project files into the container (excluding auto-07p, based on .dockerignore file)
COPY . /project

# Set working directory to /project
WORKDIR /project


