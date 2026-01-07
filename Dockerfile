# syntax=docker/dockerfile:1

FROM --platform=linux/amd64 python:3.10-bookworm

# Prevents Python from writing pyc files.
ENV PYTHONDONTWRITEBYTECODE=1
# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering.
ENV PYTHONUNBUFFERED=1

WORKDIR /app

COPY requirements requirements
RUN pip install --no-cache-dir  -r requirements/requirements.txt
RUN pip install --no-cache-dir -r requirements/requirements_gui.txt

COPY directlfq directlfq
COPY MANIFEST.in MANIFEST.in
COPY LICENSE LICENSE
COPY README.md README.md
COPY pyproject.toml pyproject.toml

RUN pip install --no-cache-dir ".[stable,gui-stable]"

ENV PORT=5006
EXPOSE 5006

# to allow other host ports than 5006
ENV BOKEH_ALLOW_WS_ORIGIN=localhost

CMD ["/usr/local/bin/directlfq", "gui", "--port", "5006"]

# build & run:
# docker build --progress=plain -t directlfq .
# DATA_FOLDER=/path/to/local/data
# docker run -p 41215:41215 -v $DATA_FOLDER:/app/data/ -t directlfq