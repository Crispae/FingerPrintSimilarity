# Use an official Python runtime as a base image
FROM python:3.9-slim

# Set environment variables to avoid buffering and to use UTF-8 encoding
ENV PYTHONUNBUFFERED=1
ENV PYTHONIOENCODING=UTF-8

# Install Java and necessary dependencies for Python
RUN apt-get update && \
    apt-get install -y default-jre && \
    rm -rf /var/lib/apt/lists/*

# Set JAVA_HOME environment variable
ENV JAVA_HOME=/usr/lib/jvm/default-java
ENV PATH="${JAVA_HOME}/bin:${PATH}"


# Set a working directory
WORKDIR /app

# Copy the entire project into the container
COPY . /app/

# Install dependencies
RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

EXPOSE 8501

CMD ["streamlit", "run", "streamlit.py", "--server.port=8501", "--server.address=0.0.0.0"]
