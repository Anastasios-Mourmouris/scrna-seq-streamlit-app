import streamlit as st
import scanpy as sc
import tempfile
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
sys.path.append("C:/Users/tasos/Desktop/TL_Project_Codes/TL_2025_scRNA-seq_pipeline_python-main")

from adata_preprocessor import adata_preprocessor, save_adata

st.set_page_config(page_title="scRNA-seq Pipeline App", layout="wide")

st.title("🔬 scRNA-seq Data Analysis Application")
st.markdown("""
Αυτή η εφαρμογή επιτρέπει την ανάλυση δεδομένων single-cell RNA-seq μέσω προκαθορισμένων βημάτων pre-processing και visualizations.
""")

# Sidebar for navigation
st.sidebar.title("Μενού")
section = st.sidebar.radio("Πλοήγηση",
                           ["Αρχική", "Upload Δεδομένων", "Preprocessing", "Οπτικοποιήσεις", "Διαφορική Έκφραση (DEG)",
                            "Ομάδα"])

# Home Tab
if section == "Αρχική":
    st.header("📘 Καλώς ήρθατε!")
    st.markdown("""
    Η εφαρμογή αυτή υποστηρίζει την ανάλυση δεδομένων μονοκυτταρικής έκφρασης RNA (scRNA-seq).

    **Βήματα υποστηριζόμενα:**
    - Μεταφόρτωση αρχείου .h5ad
    - Καθαρισμός δεδομένων (π.χ. φιλτράρισμα γονιδίων και κυττάρων)
    - Προετοιμασία για DEG και οπτικοποιήσεις
    - Οπτικοποίηση γονιδιακής έκφρασης
    - Ανάλυση διαφορικής έκφρασης

    """)

# Upload Tab
elif section == "Upload Δεδομένων":
    st.header("📤 Μεταφόρτωση Αρχείου .h5ad")
    uploaded_file = st.file_uploader("Επέλεξε αρχείο .h5ad", type=["h5ad"])

    if uploaded_file:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
            tmp.write(uploaded_file.read())
            st.session_state["adata_path"] = tmp.name
            st.success("Το αρχείο ανέβηκε επιτυχώς και αποθηκεύτηκε προσωρινά.")

# Preprocessing Tab
elif section == "Preprocessing":
    st.header("⚙️ Προεπεξεργασία Δεδομένων")

    if "adata_path" not in st.session_state:
        st.warning("Παρακαλώ πρώτα ανέβασε αρχείο .h5ad στην αντίστοιχη καρτέλα.")
    else:
        st.markdown("Ορισμός παραμέτρων για preprocessing:")

        min_genes = st.slider("min_genes", 0, 500, 100)
        min_cells = st.slider("min_cells", 1, 10, 3)
        n_genes_min = st.slider("n_genes_min", 500, 3000, 1000)
        n_genes_max = st.slider("n_genes_max", 5000, 20000, 10000)
        n_counts_max = st.slider("n_counts_max", 10000, 100000, 30000)
        pc_mito = st.slider("% Μιτοχονδριακών", 0, 100, 20)
        pc_rib = st.slider("% Ριβοσωμικών", 0, 100, 25)

        if st.button("🔄 Εκτέλεση Preprocessing"):
            adata_filtered = adata_preprocessor(
                st.session_state["adata_path"],
                min_genes=min_genes,
                min_cells=min_cells,
                n_genes_min=n_genes_min,
                n_genes_max=n_genes_max,
                n_counts_max=n_counts_max,
                pc_mito=pc_mito,
                pc_rib=pc_rib,
                debug=True
            )

            if adata_filtered.shape[0] == 0 or adata_filtered.shape[1] == 0:
                st.error("❌ Δεν υπάρχουν δεδομένα μετά το preprocessing. Δοκίμασε χαλαρότερα φίλτρα.")
                st.stop()

            output_path = os.path.join(tempfile.gettempdir(), "filtered_adata.h5ad")
            save_adata(adata_filtered, output_path)
            st.session_state["adata_filtered"] = adata_filtered
            st.session_state["output_path"] = output_path

            st.success("Preprocessing ολοκληρώθηκε. Τα δεδομένα φορτώθηκαν στη μνήμη και αποθηκεύτηκαν προσωρινά.")

            with open(output_path, "rb") as file:
                st.download_button(label="⬇️ Κατέβασε το φιλτραρισμένο .h5ad",
                                   data=file,
                                   file_name="filtered_adata.h5ad")



# Visualization Tab
elif section == "Οπτικοποιήσεις":
    st.header("📊 Οπτικοποιήσεις Έκφρασης Γονιδίων")

    if "adata_filtered" not in st.session_state:
        st.warning("Πρέπει να εκτελέσεις πρώτα το preprocessing.")
    else:
        adata = st.session_state["adata_filtered"]

        if adata.shape[0] == 0:
            st.error("❌ Δεν υπάρχουν κύτταρα για plotting UMAP. Δοκίμασε διαφορετικές παραμέτρους preprocessing.")
            st.stop()

        st.subheader("Διαδραστικό UMAP")
        if "X_umap" not in adata.obsm:
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)

        fig_umap = sc.pl.umap(adata, color=["n_genes", "n_counts"], return_fig=True, show=False)
        st.pyplot(fig_umap)

        st.subheader("Επιλογή Γονιδίων για Plot")
        gene_list = st.text_input("Δώσε ονόματα γονιδίων χωρισμένα με κόμμα:", value="CD3D,MS4A1")
        genes = [g.strip() for g in gene_list.split(",")]

        if st.button("📈 Plot Gene Expression"):
            try:
                import matplotlib.pyplot as plt

                sc.pl.violin(
                    adata,
                    keys=genes,
                    groupby="leiden" if "leiden" in adata.obs else None,
                    jitter=0.4,
                    multi_panel=True,
                    show=False  # ΜΗΝ εμφανίσεις με plt.show()
                )

                st.pyplot(plt.gcf())  # ⬅️ Πάρε το ενεργό figure και εμφάνισέ το στο Streamlit

            except Exception as e:
                st.error(f"Σφάλμα στην οπτικοποίηση: {e}")

# DEG Tab
elif section == "Διαφορική Έκφραση (DEG)":
    st.header("📈 Ανάλυση Διαφορικής Έκφρασης")

    if "adata_filtered" not in st.session_state:
        st.warning("Απαιτείται φιλτραρισμένο αρχείο .h5ad από preprocessing.")
    else:
        adata = st.session_state["adata_filtered"]

        if adata.shape[0] == 0:
            st.error("❌ Δεν υπάρχουν κύτταρα στο adata για DEG.")
            st.stop()

        if "leiden" not in adata.obs:
            st.info("Εκτελείται clustering μέσω Leiden...")
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.leiden(adata)
            st.success("Έγινε clusterization με βάση το Leiden.")

        groupby = "leiden"

        st.markdown("Παραγωγή διαφορικής έκφρασης γονιδίων ανά cluster:")
        if st.button("🚀 Υπολογισμός DEG"):
            try:
                sc.tl.rank_genes_groups(adata, groupby=groupby, method='t-test')
                result = adata.uns['rank_genes_groups']
                groups = result['names'].dtype.names

                deg_table = pd.DataFrame({
                    group + '_' + key: result[key][group]
                    for group in groups
                    for key in ['names', 'scores', 'logfoldchanges', 'pvals_adj']
                })

                st.session_state["deg_table"] = deg_table
                st.dataframe(deg_table.head(20))

                csv = deg_table.to_csv(index=False).encode('utf-8')
                st.download_button("📥 Κατέβασε DEG .csv", csv, "deg_results.csv", "text/csv")

                # Volcano Plot (χωρίς return_fig)
                with st.expander("📊 Προβολή Volcano Plot (top 25 γονίδια)"):
                    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False)
                    st.pyplot(plt.gcf())
                    plt.clf()

            except Exception as e:
                st.error(f"❌ Σφάλμα κατά τον υπολογισμό ή την εμφάνιση των DEG: {e}")

# Team Tab
elif section == "Ομάδα":
    st.header("👨‍💻 Πληροφορίες Ομάδας")
    st.markdown("""
    - **Όνομα Φοιτητή/Ομάδας:** 
    - **Μέλη:**
      - Αναστάσιος Μουρμούρης 
       - **Συνεισφορά:**
            - Ολοκλήρωση της ανάλυσης DEG με Scanpy\n
            - Προβολή αποτελεσμάτων: πίνακες, volcano\n
            - Υλοποίηση επιλογής γονιδίων και violin plots\n
            - Χειρισμός session state για αποθήκευση αποτελεσμάτων\n
      - Στέφανος Λάμπρου
       - **Συνεισφορά:**
            - Υλοποίηση adata_preprocessor.py\n
            - Ανάγνωση και αποθήκευση .h5ad\n
            - Προσωρινή διαχείριση αρχείων με tempfile\n
            - Σύνδεση με Streamlit tabs για preprocessing\n
            - Έλεγχος για άδεια δεδομένα / σφάλματα\n
      - Ορφέας Λάμπρου
       - **Συνεισφορά:**\n
             - Σχεδιασμός tab-based Streamlit UI (main.py)\n
             - Οπτική διάταξη στοιχείων (sliders, inputs, buttons)\n
             - Δημιουργία Dockerfile, requirements.txt, τεκμηρίωση χρήσης\n
             - UML διαγράμματα (Use Case, Class)\n
             - Συγγραφή μεγάλου μέρους του τελικού LaTeX report\n
    """)
