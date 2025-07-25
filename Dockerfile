FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Copy the rest of the application files
COPY . .

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libblas-dev \
    liblapack-dev \
    libgfortran5 \
    git \
    && rm -rf /var/lib/apt/lists/*


# Install dependencies
RUN pip install -e .

# Expose default Node.js port
EXPOSE 8000

# Start the application
CMD ["python", "main.py"]