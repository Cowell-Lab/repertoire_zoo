from connexion import AsyncApp
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

import repertoire_zoo.bird as br
import repertoire_zoo.giraffe as gr

async def post_greeting(name: str):
    return f"Hello {name}\n", 200

async def gene_usage(body):
    try:
        # Validate input
        data_path = body["data_path"]
        repertoire_id = body["repertoire_id"]
        processing_stage = body["processing_stage"]
    except (TypeError, KeyError) as e:
        return {"error": f"Missing or invalid input: {e}"}, 400

    try:
        # Load data
        gene_usage_df = br.load_gene_usage_data(
            repcalc_dir=data_path,
            repertoire_id=repertoire_id,
            processing_stage=processing_stage
        )
        print(gene_usage_df.columns)
        print(gene_usage_df.head())
        
    except Exception as e:
        return {"error": f"Failed to load gene usage data: {e}"}, 500

    try:
        # Generate plot
        fig, ax = gr.plot_duplicate_frequency(gene_usage_df)
        # Save plot
        filename = f"gene_usage_{repertoire_id}.png"
        file_path = f"/Figures/{filename}"
        fig.savefig(file_path, format='png')
        print(f"Figure save in {file_path}")
        plt.close(fig)
        return {"image_path": str(file_path)}, 200, {"Content-Type": "application/json"}

    except Exception as e:
        return {"error": f"Failed to generate or save plot: {e}"}, 500

# Set up the Connexion app
app = AsyncApp(__name__)
app.add_api("/repertoire_zoo/openapi/vdjserver-image.yaml")