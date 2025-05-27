# scRNAseqApp

A Streamlit-based interactive application for single-cell RNA-seq data analysis.

## 🔧 Requirements

- Docker installed
- Input file in `.h5ad` format

## 🚀 Usage

```bash
docker build -t scrna-seq-app .
docker run -p 8501:8501 scrna-seq-app
```

Then open your browser at [http://localhost:8501](http://localhost:8501)

## 📁 Files

- `main.py`: Main Streamlit interface
- `adata_preprocessor.py`: Preprocessing pipeline
- `Dockerfile`: Container definition
- `requirements.txt`: Python dependencies