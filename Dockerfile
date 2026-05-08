FROM rocker/r-ver:4.4.2

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    make \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
  && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("jsonlite", "testthat", "ggplot2", "patchwork", "reshape2", "scales", "pROC"), repos = "https://cloud.r-project.org")'

COPY requirements-carf.txt /tmp/requirements-carf.txt
RUN python3 -m venv /opt/carf-venv
ENV PATH="/opt/carf-venv/bin:${PATH}"
RUN pip install --no-cache-dir -r /tmp/requirements-carf.txt

WORKDIR /workspace
COPY . /workspace

CMD ["make", "carf-v1"]
