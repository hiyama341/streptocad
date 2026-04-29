# To run a build that works on AWS us the following command: 
#docker build --platform linux/amd64 -t streptocad .

# Use Python 3.11 as the base image
FROM python:3.11-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    libc6-dev \
    libffi-dev \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libsqlite3-dev \
    libgdbm-dev \
    libdb5.3-dev \
    liblzma-dev \
    libexpat1-dev \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /web_app

# Copy dependency metadata first to make Docker layer caching useful.
COPY pyproject.toml README.md LICENSE ./
COPY streptocad ./streptocad

# Upgrade pip, setuptools, and wheel.
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

RUN pip install --no-cache-dir ".[app]"

# Copy the rest of the application code into the container
COPY . .

# Install the current project
RUN pip install --no-cache-dir -e ".[app]"

# Expose the port that Dash runs on (default is 8050)
EXPOSE 8050

# Set environment variables
ENV FLASK_ENV=production
ENV DOCKER_CONTAINER=true

# Command to run the application
CMD ["python", "web_app/application.py"]
