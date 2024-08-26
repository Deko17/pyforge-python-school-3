FROM continuumio/miniconda3
RUN conda install -c conda-forge rdkit -y

WORKDIR /pyforge-python-school-3

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY src/ src/

CMD ["uvicorn", "src.main:app", "--host", "0.0.0.0", "--port", "8000"]
