#/Users/johannawahn/Desktop/iGEM/Modeling
# -*- coding: utf-8 -*-
"""
Created on %(29/6/2020)s

@author: %(johananWahn)s
"""
from dash import Dash, html, dcc, callback, Output, Input, State
import plotly.graph_objects as go
import pandas as pd
import random
import plotly.express as px



data = pd.read_csv('coordinates_GroupA.csv')
data = data.rename(columns={'Precursor.Id':'peptide'})

samples = data.iloc[:,1].unique().tolist()

DEFAULT_SIZE = 4
HIGHLIGHT_SIZE = 10
SELECTED = None

get_colors = lambda n: ["#%06x" % random.randint(0, 0xFFFFFF) for _ in range(n)]
palette = get_colors(len(data['ProteinGroup'].unique()))

styles = {}
for nbr, unique_prot in enumerate(data['ProteinGroup'].unique()):
    styles[unique_prot] = dict(marker=dict(color=palette[nbr], size=DEFAULT_SIZE))

app = Dash(__name__)

app.layout = html.Div([
    html.H1(children='Visualization of proteomic run', style={'textAlign':'center'}),
    dcc.Dropdown(samples, samples[2] , id='dropdown-selection'),
    dcc.Graph(id='graph-content')
])


@app.callback(
    Output('graph-content', 'figure'),
    Input('graph-content', 'hoverData'), 
    Input('dropdown-selection', 'value'),
    State('graph-content', 'figure'),
)

def update_graph(hover_data, value, fig):
    
    df = data[data['Run']==value]
    
    global SELECTED
    if hover_data is not None:
        if SELECTED is not None:
            fig['data'][SELECTED]['marker']['size'] = DEFAULT_SIZE

        ProteinGroup = hover_data['points'][0]['customdata']
        SELECTED = ix = [scatter['name'] for scatter in fig['data']].index(ProteinGroup)
        fig['data'][ix]['marker']['size'] = HIGHLIGHT_SIZE
        
    else:
        fig = go.Figure()
        fig.update_xaxes(showgrid=False, range=(0, 2250))
        fig.update_yaxes(showgrid=False, range=(0, 2287), 
                         scaleanchor="x",scaleratio=1,)
        fig.update_layout(autosize=True,
                          xaxis_title = 'RT (sec)',
                          yaxis_title = 'mz',
                          legend_title = 'Protein Groups',
                          title = value) 


        
        for ProteinGroup, peptides in df.groupby('ProteinGroup'):
            fig.add_scatter(
                x=peptides.RT, #RT
                y=peptides.mz, #mz
                name=ProteinGroup, #protein
                customdata=peptides.ProteinGroup, #protein for the peptide
                mode='markers',
                hovertext=peptides.peptide,
                **styles[ProteinGroup]
                
            )


        fig.update_layout(height=1000, hovermode='closest', uirevision='static')
    return fig




if __name__ == '__main__':
    app.run_server(debug=True)
