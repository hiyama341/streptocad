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

# Copy the requirements file into the container
COPY requirements.txt .

# Upgrade pip, setuptools, and wheel
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Install the required packages, excluding the -e . line
RUN grep -v '^-e' requirements.txt > requirements_filtered.txt && \
    pip install --no-cache-dir -r requirements_filtered.txt

# Copy the rest of the application code into the container
COPY . .

# Install the current project
RUN pip install -e .

# Expose the port that Dash runs on (default is 8050)
EXPOSE 8050

# Set environment variables
ENV FLASK_ENV=production
ENV DOCKER_CONTAINER=true

# Command to run the application
CMD ["python", "web_app/application.py"]
