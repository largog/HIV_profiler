import dash
from dash import dcc, html, dash_table, Input, Output, State
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import flask

# Importar pipeline
from pipeline_hiv import full_pipeline

app = dash.Dash(__name__)
app.title = "Dashboard VIH – Farmacología"
server = app.server

# Servir imágenes desde /results
@server.route("/results/<path:path>")
def serve_results(path):
    return flask.send_from_directory("results", path)


# ==========================================================
#                     LAYOUT
# ==========================================================

app.layout = html.Div(
    style={"fontFamily": "Arial", "margin": "20px"},
    children=[
        html.H1("Análisis de variantes de VIH", style={"textAlign": "center"}),

        html.Hr(),

        # ----------------- ENTRADA -----------------
        html.H3("1. Entrada de datos"),
        html.Div(
            style={"display": "flex", "gap": "10px", "flexWrap": "wrap"},
            children=[
                html.Label("Ruta archivo FASTA:"),
                dcc.Input(
                    id="query-path",
                    type="text",
                    value="data_test/B1.fasta",
                    style={"width": "35%"},
                ),
                html.Label("Base de datos BLAST:"),
                dcc.Input(
                    id="db-path",
                    type="text",
                    value="database/HIV1_REF",
                    style={"width": "35%"},
                ),
                html.Button("Ejecutar análisis", id="run-button", n_clicks=0),
            ],
        ),

        html.Hr(),

        # ----------------- BLAST -----------------
        html.H3("2. Resultado del BLAST y subtipo detectado"),

        html.Div(id="summary-cards", style={"display": "flex", "gap": "20px"}),

        html.Br(),

        html.Div(
            style={"display": "flex", "gap": "20px"},
            children=[
                dcc.Graph(id="identity-graph", style={"width": "50%"}),
                dcc.Graph(id="alignment-graph", style={"width": "50%"}),
            ],
        ),

        html.Br(),

        html.H4("Tabla de resultado del BLAST"),
        dash_table.DataTable(
            id="blast-table",
            columns=[],
            data=[],
            row_selectable="single",
            style_cell={"textAlign": "center"},
            style_header={"fontWeight": "bold"},
        ),

        html.Hr(),

        # ----------------- PHARMA -----------------
        html.H3("3. Datos farmacológicos del subtipo"),
        dash_table.DataTable(
            id="pharma-table",
            columns=[],
            data=[],
            page_size=5,
            style_header={"fontWeight": "bold"},
            style_cell={"textAlign": "center"},
        ),

        html.Br(),

        # ----------------- TOP 3 -----------------
        html.H3("Top 3 fármacos más potentes"),
        dash_table.DataTable(
            id="top3-table",
            columns=[],
            data=[],
            row_selectable="single",
            style_header={"fontWeight": "bold"},
            style_cell={"textAlign": "center"},
        ),

        html.Br(),

        # ----------------- MOLÉCULA -----------------
        html.H3("Visualización del fármaco seleccionado"),

        html.Div(
            style={"display": "flex", "gap": "20px"},
            children=[
                html.Img(
                    id="mol-image",
                    style={"width": "300px", "border": "1px solid #ccc"},
                ),
                dash_table.DataTable(
                    id="mol-properties",
                    columns=[],
                    data=[],
                    style_cell={"textAlign": "left"},
                    style_header={"fontWeight": "bold"},
                ),
            ],
        ),

        # Stores (memoria)
        dcc.Store(id="store-pharma"),
        dcc.Store(id="store-top3"),
        dcc.Store(id="store-blast"),
    ],
)


# ==========================================================
#                CALLBACK 1 — EJECUTAR PIPELINE
# ==========================================================

@app.callback(
    [
        Output("summary-cards", "children"),
        Output("identity-graph", "figure"),
        Output("alignment-graph", "figure"),
        Output("blast-table", "data"),
        Output("blast-table", "columns"),
        Output("store-pharma", "data"),
        Output("store-top3", "data"),   # ← ahora contiene rdkit_info
        Output("store-blast", "data"),
    ],
    Input("run-button", "n_clicks"),
    State("query-path", "value"),
    State("db-path", "value"),
    prevent_initial_call=True,
)
def run_analysis(n_clicks, query_path, db_path):

    blast_info, pharma_df, top3, rdkit_info = full_pipeline(
        query_path, db_path, data_folder="."
    )

    # ---------- TARJETAS ----------
    cards = [
        html.Div(
            style={"border": "1px solid #ddd", "padding": "10px"},
            children=[html.H4("Subtipo"), html.P(blast_info["subtype"])],
        ),
        html.Div(
            style={"border": "1px solid #ddd", "padding": "10px"},
            children=[html.H4("Identidad"), html.P(f"{blast_info['percent_identity']}%")],
        ),
        html.Div(
            style={"border": "1px solid #ddd", "padding": "10px"},
            children=[html.H4("Alineamiento"), html.P(blast_info["alignment_length"])],
        ),
        html.Div(
            style={"border": "1px solid #ddd", "padding": "10px"},
            children=[html.H4("Cobertura"),html.P(f"{blast_info['coverage']}")],
        )
    ]

    # ---------- GRAFICA IDENTIDAD y COBERTURA ----------
    coverage_p = blast_info["coverage"] * 100

    fig_identity = px.bar(
        x=["Identidad (%)", "Cobertura (%)"],
        y=[blast_info["percent_identity"], coverage_p],
        text=[blast_info["percent_identity"], coverage_p],
        range_y=[0, 100],
        title="Porcentaje de identidad y cobertura",
    )

    # ---------- GRAFICA ALIGNMENT ESTILO BLAST ----------
    alignment_df = pd.DataFrame([
        {"Secuencia": "Referencia", "Inicio": 0, "Fin": blast_info["reference_length"]},
        {"Secuencia": "Alineamiento", "Inicio": blast_info["S_start"], "Fin": blast_info["S_end"]},
    ])

    fig_align = go.Figure()

    # Barra de referencia (0 → reference_length)
    fig_align.add_trace(go.Bar(
        y=["Referencia"],
        x=[blast_info["reference_length"]],
        orientation="h",
        name="Referencia",
        base=[0],
        text=[blast_info["reference_length"]],
    ))

    # Barra de alineamiento (S_start → S_end)
    fig_align.add_trace(go.Bar(
        y=["Alineamiento"],
        x=[blast_info["S_end"] - blast_info["S_start"]],
        orientation="h",
        name="Alineamiento",
        base=[blast_info["S_start"]],
        text=[blast_info["S_end"]],
    ))

    fig_align.update_layout(
        title="Alineamiento",
        barmode="overlay",
        xaxis_title="Posición en el genoma",
        yaxis_title="Secuencia",
        showlegend = False
    )


    # ---------- TABLA BLAST ----------
    blast_df = pd.DataFrame([blast_info])

    return (
        cards,
        fig_identity,
        fig_align,
        blast_df.to_dict("records"),
        [{"name": c, "id": c} for c in blast_df.columns],
        pharma_df.to_dict("records"),
        rdkit_info.to_dict("records"),   # ← SOLO rdkit_info tiene “Image”
        blast_df.to_dict("records")[0],
    )


# ==========================================================
#        CALLBACK 2 — MOSTRAR PHARMA + TOP3
# ==========================================================

@app.callback(
    [
        Output("pharma-table", "data"),
        Output("pharma-table", "columns"),
        Output("top3-table", "data"),
        Output("top3-table", "columns"),
    ],
    Input("blast-table", "active_cell"),
    State("store-pharma", "data"),
    State("store-top3", "data"),
    prevent_initial_call=True,
)
def show_pharma_top3(active_cell, pharma_data, top3_data):

    if active_cell is None:
        return [], [], [], []

    # ------- FILTRAR columnas farmacológicas -------
    pharma_df = pd.DataFrame(pharma_data)
    cols_to_show = ["BindingDB Reactant_set_id", "Ligand SMILES", "Ki (nM)"]
    pharma_df = pharma_df[cols_to_show]

    # ------- CONVERTIR top3_data A DATAFRAME -------
    top3_df = pd.DataFrame(top3_data)

    # ------- TOP 3 MOSTRAR RANKING + NOMBRE --------
    ranking_df = pd.DataFrame({
        "Ranking": [1, 2, 3],
        "Fármaco": top3_df["Commercial Name"].tolist()
    })

    return (
        pharma_df.to_dict("records"),
        [{"name": c, "id": c} for c in cols_to_show],

        ranking_df.to_dict("records"),
        [
            {"name": "Ranking", "id": "Ranking"},
            {"name": "Fármaco", "id": "Fármaco"}
        ],
    )


# ==========================================================
#        CALLBACK 3 — IMAGEN + PROPIEDADES MOLÉCULA
# ==========================================================

@app.callback(
    [
        Output("mol-image", "src"),
        Output("mol-properties", "data"),
        Output("mol-properties", "columns"),
    ],
    Input("top3-table", "active_cell"),
    State("store-top3", "data"),
    prevent_initial_call=True,
)
def display_molecule(active_cell, top3_data):

    if active_cell is None or top3_data is None:
        return None, [], []

    idx = active_cell["row"]
    mol_info = top3_data[idx]

    # Imagen RDKit
    src = f"/{mol_info['Image']}"

    # Tabla de propiedades
    prop_dict = {
        "Commercial Name": mol_info.get("Commercial Name"),
        "IUPAC Name": mol_info.get("IUPAC Name"),
    }

    # luego agrega las demás propiedades
    for k, v in mol_info.items():
        if k not in ["Image", "SMILES", "Commercial Name", "IUPAC Name"]:
            prop_dict[k] = v

    prop_table = [{"Propiedad": k, "Valor": v} for k, v in prop_dict.items()]
    columns = [{"name": "Propiedad", "id": "Propiedad"},
               {"name": "Valor", "id": "Valor"}]

    return src, prop_table, columns


# ==========================================================
#                    RUN SERVER
# ==========================================================

if __name__ == "__main__":
    app.run(debug=True)

