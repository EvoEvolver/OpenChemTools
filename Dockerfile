FROM ghcr.io/astral-sh/uv:python3.12-alpine

# Set working directory
WORKDIR /app

# Copy the rest of the application files
COPY . .

RUN apk add --no-cache \
    git \
    cmake \
    make \
    gcc \
    g++ \
    gfortran \
    pkgconfig \
    hdf5-dev \
    openblas-dev \
    lapack-dev \
    fftw-dev \
    musl-dev


# Install dependencies
RUN uv pip install --system -e .
RUN uv pip install --system --prefer-binary pyscf

# Expose default Node.js port
EXPOSE 8000

# Start the application
CMD ["python", "main.py"]