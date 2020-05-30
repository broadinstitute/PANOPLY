FROM broadgdac/run-r:3.3.2

# Install any libraries necessary to run the task
# RUN set -ex \
#     && apt-get update \
#     && apt-get install -y --no-install-recommends \
#         libcairo2-dev \
#         libxt-dev \
#     && rm -rf /var/lib/apt/lists/* \
#     && install2.r -e \
#         Cairo \
#         data.table \
#     && rm -rf /tmp/*

# Copy the built tool and any supporting files into the image
COPY src /src

# Set the working directory
WORKDIR src
