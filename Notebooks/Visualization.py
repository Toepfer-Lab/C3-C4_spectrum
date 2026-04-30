import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import sys
sys.path.append("../Functions")
from Auxiliary import *
from Flux import *
import numpy as np
import altair as alt
from cobra.io import read_sbml_model
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
pd.options.mode.chained_assignment = None  # default='warn'


def comparison_ATP_barplots_3D(df_C3, df_C4, df_C4_only, PPFD, Vc_Vo, CO2, PC_BS, NH3):
    """
    Possible values for selected parameters (initial parameter scans):

    PPFD = [100, 200, 300, 400, 500, 600, 700, 800, 1000]
    CO2 = [20, 1000]
    Vc_Vo = [50, 55, 60 ,65, 70, 75, 80, 85, 90, 95, 99, 100]
    PC_BS = [1, .1]

    """

    selected_C3 = select_data(prepared_df=df_C3, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3=NH3)
    selected_C4 = select_data(prepared_df=df_C4, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3=NH3)
    selected_C4_only = select_data(prepared_df=df_C4_only, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3=NH3)

    data_C4 = selected_C4.copy()
    data_C4.index=data_C4.V_BS_M
    data_C4["Classification"] = "Variable"

    data_C3 = selected_C3.copy()
    data_C3.index=data_C3.V_BS_M
    data_C3["Classification"] = "Only C3"

    data_C4_only = selected_C4_only.copy()
    data_C4_only.index=data_C4_only.V_BS_M
    data_C4_only["Classification"] = "Only C4"

    df = pd.concat([data_C4, data_C3, data_C4_only])

    symbol_mapping = {
        'C3': 'circle',
        'Unidentified': 'cross',
        'C3_C4_Type_I': 'square',
        'C3_C4_Type_II': 'diamond',
        'C4': 'triangle-up'
    }

    # Create a new column in the DataFrame that maps Photosynthesis to symbols
    df['symbol'] = df['Photosynthesis'].map(symbol_mapping)

    
    # Creating an adjusted ATP value for the point chart
    df['Adjusted_ATP'] = df['Total_ATP'] + (df['Total_ATP'].max() * 0.05)

    # Creating the bar chart
    bar = alt.Chart(df).mark_bar().encode(
        x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
        y=alt.Y('Total_ATP:Q', axis=alt.Axis(grid=False, title="ATP Turnover (µmol m-2 s-1)")),
        color=alt.Color('Classification:N', scale=alt.Scale(range=['#1a1f16', '#1E3F20', "#4A7856"]))
    )

    # Creating the point chart with a custom shape legend
    symbols = alt.Chart(df).mark_point(size=100).encode(
        x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
        y=alt.Y('Adjusted_ATP:Q', axis=alt.Axis(grid=False)),
        shape=alt.Shape(
            'Photosynthesis:N',  # Use the Photosynthesis column to map to the shapes
            scale=alt.Scale(
                domain=['C3', 'Unidentified', 'C3_C4_Type_I', 'C3_C4_Type_II', 'C4'],
                range=['circle', 'cross', 'square', 'diamond', 'triangle-up']
            ),
            legend=alt.Legend(title="Photosynthesis Type")  # Custom legend title
        ),
        color=alt.Color('Classification:N', scale=alt.Scale(range=['#C1A32A', '#28B34A']))
    )

    # Combining the bar and point charts using facet to ensure a shared y-axis
    chart = alt.layer(bar, symbols).facet(
        column=alt.Column('V_BS_M:N', title="Volume Ratio")
    ).resolve_scale(
        y='shared'  
    ).configure_view(
        strokeOpacity=0  # Remove grid lines around column clusters
    )

    return chart

def comparison_ATP_split_barplots_3D(df_C3, df_C4, df_C4_only, PPFD, Vc_Vo, CO2, PC_BS, NH3, min_y, max_y):

    volume = [0.1, 0.5, 1, 1.5, 2]
    
    selected_C3 = select_data(prepared_df=df_C3, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
    selected_C4 = select_data(prepared_df=df_C4, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
    selected_C4_only = select_data(prepared_df=df_C4_only, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)

    selected_C3 = selected_C3[selected_C3.V_BS_M.isin(volume)]
    selected_C4 = selected_C4[selected_C4.V_BS_M.isin(volume)]
    selected_C4_only = selected_C4_only[selected_C4_only.V_BS_M.isin(volume)]


    # Use melt to reshape the data
    melted_df_C4 = pd.melt(
        selected_C4,
        id_vars=['V_BS_M', 'Photosynthesis'],
        value_vars=['ATP M', 'ATP BS'],
        var_name='ATP',
        value_name='Value'
    )

    melted_df_C3 = pd.melt(
        selected_C3,
        id_vars=['V_BS_M', 'Photosynthesis'],
        value_vars=['ATP M', 'ATP BS'],
        var_name='ATP',
        value_name='Value'
    )

    melted_df_C4_only = pd.melt(
        selected_C4_only,
        id_vars=['V_BS_M', 'Photosynthesis'],
        value_vars=['ATP M', 'ATP BS'],
        var_name='ATP',
        value_name='Value'
    )

    melted_df_C3["Classification"] = "Only C3"
    melted_df_C4["Classification"] = "Variable"
    melted_df_C4_only["Classification"] = "Only C4"

    df = pd.concat([melted_df_C3, melted_df_C4, melted_df_C4_only])

    symbol_mapping = {
        'C3': 'circle',
        'Unidentified': 'cross',
        'C3_C4_Type_I': 'square',
        'C3_C4_Type_II': 'diamond',
        'C4': 'triangle-up'
    }

    # Create a new column in the DataFrame that maps Photosynthesis to symbols
    df['symbol'] = df['Photosynthesis'].map(symbol_mapping)


    # Creating an adjusted ATP value for the point chart
    df['Adjusted_Value'] = df['Value'] + (df['Value'].max() * 0.03)

    # Creating the bar chart
    bar = alt.Chart(df).mark_bar().encode(
        x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
        y=alt.Y('sum(Value):Q', axis=alt.Axis(grid=False, title="ATP Turnover (µmol m-2 s-1)"),
               scale=alt.Scale(domain=[min_y, max_y])
               ),
        color=alt.Color('ATP:N', scale=alt.Scale(range=['#1a1f16', '#4A7856']))
    )

    # Creating the point chart with a custom shape legend
    symbols = alt.Chart(df).mark_point(size=100).encode(
        x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
        y=alt.Y('sum(Adjusted_Value):Q', axis=alt.Axis(grid=False),
               scale=alt.Scale(domain=[min_y, max_y])
               ),
        shape=alt.Shape(
            'Photosynthesis:N',  # Use the Photosynthesis column to map to the shapes
            scale=alt.Scale(
                domain=['C3', 'Unidentified', 'C3_C4_Type_I', 'C3_C4_Type_II', 'C4'],
                range=['circle', 'cross', 'square', 'diamond', 'triangle-up']
            ),
            legend=alt.Legend(title="Photosynthesis Type")  
        )
    )

    # Combining the bar and point charts using facet to ensure a shared y-axis
    chart = alt.layer(bar, symbols).facet(
        column=alt.Column('V_BS_M:N', title="Volume Ratio")
    ).resolve_scale(
        y='shared'  
    ).configure_view(
        strokeOpacity=0  # Remove grid lines around column clusters
    )

    return chart

def comparison_Biomass_barplots_3D(df_C3, df_C4, df_C4_only, PPFD, Vc_Vo, CO2, PC_BS, NH3):
    """
    Possible values for selected parameters (initial parameter scans):

    PPFD = [100, 200, 300, 400, 500, 600, 700, 800, 1000]
    CO2 = [20, 1000]
    Vc_Vo = [50, 55, 60 ,65, 70, 75, 80, 85, 90, 95, 99, 100]
    PC_BS = [1, .1]

    """
    volume = [0.1, 0.5, 1, 1.5, 2]
    
    selected_C3 = select_data(prepared_df=df_C3, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
    selected_C4 = select_data(prepared_df=df_C4, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
    selected_C4_only = select_data(prepared_df=df_C4_only, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)

    selected_C3 = selected_C3[selected_C3.V_BS_M.isin(volume)]
    selected_C4 = selected_C4[selected_C4.V_BS_M.isin(volume)]
    selected_C4_only = selected_C4_only[selected_C4_only.V_BS_M.isin(volume)]

    data_C4 = selected_C4.copy()
    data_C4.index=data_C4.V_BS_M
    data_C4["Classification"] = "Variable"

    data_C3 = selected_C3.copy()
    data_C3.index=data_C3.V_BS_M
    data_C3["Classification"] = "Only C3"

    data_C4_only = selected_C4_only.copy()
    data_C4_only.index=data_C4_only.V_BS_M
    data_C4_only["Classification"] = "Only C4"

    df = pd.concat([data_C4, data_C3, data_C4_only])

    symbol_mapping = {
        'C3': 'circle',
        'Unidentified': 'cross',
        'C3_C4_Type_I': 'square',
        'C3_C4_Type_II': 'diamond',
        'C4': 'triangle-up'
    }

    # Create a new column in the DataFrame that maps Photosynthesis to symbols
    df['symbol'] = df['Photosynthesis'].map(symbol_mapping)

    
    # Creating an adjusted ATP value for the point chart
    df['Adjusted_Biomass'] = df['Biomass'] + (df['Biomass'].max() * 0.05)

    # Creating the bar chart
    bar = alt.Chart(df).mark_bar().encode(
        x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
        y=alt.Y('Biomass:Q', axis=alt.Axis(grid=False, title="Biomass Production (µmol m-2 s-1)")),
        color=alt.Color('Classification:N', scale=alt.Scale(range=['#721817 ', '#FA9F42', "#2B4162"]))
    )

    # Creating the point chart with a custom shape legend
    symbols = alt.Chart(df).mark_point(size=100).encode(
        x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
        y=alt.Y('Adjusted_Biomass:Q', axis=alt.Axis(grid=False)),
        shape=alt.Shape(
            'Photosynthesis:N',  # Use the Photosynthesis column to map to the shapes
            scale=alt.Scale(
                domain=['C3', 'Unidentified', 'C3_C4_Type_I', 'C3_C4_Type_II', 'C4'],
                range=['circle', 'cross', 'square', 'diamond', 'triangle-up']
            ),
            legend=alt.Legend(title="Photosynthesis Type")  # Custom legend title
        ),
        color=alt.Color('Classification:N', scale=alt.Scale(range=['#C03924  ', '#1F7A99 ']))
    )

    # Combining the bar and point charts using facet to ensure a shared y-axis
    chart = alt.layer(bar, symbols).facet(
        column=alt.Column('V_BS_M:N', title="Volume Ratio")
    ).resolve_scale(
        y='shared'  # Ensures that the y-axis is the same across all columns
    ).configure_view(
        strokeOpacity=0  # Remove grid lines around column clusters
    )

    return chart

def comparison_Light_barplots_3D(df_C3, df_C4, df_C4_only, PPFD, Vc_Vo, CO2, PC_BS, NH3, min_y, max_y):

        volume = [0.1, 0.5, 1, 1.5, 2]
        
        selected_C3 = select_data(prepared_df=df_C3, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
        selected_C4 = select_data(prepared_df=df_C4, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
        selected_C4_only = select_data(prepared_df=df_C4_only, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)

        selected_C3 = selected_C3[selected_C3.V_BS_M.isin(volume)]
        selected_C4 = selected_C4[selected_C4.V_BS_M.isin(volume)]
        selected_C4_only = selected_C4_only[selected_C4_only.V_BS_M.isin(volume)]

        # Use melt to reshape the data
        melted_df_C4 = pd.melt(
            selected_C4,
            id_vars=['V_BS_M', 'Photosynthesis'],
            value_vars=['Light Mesophyll', 'Light Bundle Sheath'],
            var_name='Light_Cell',
            value_name='Value'
        )

        melted_df_C3 = pd.melt(
            selected_C3,
            id_vars=['V_BS_M', 'Photosynthesis'],
            value_vars=['Light Mesophyll', 'Light Bundle Sheath'],
            var_name='Light_Cell',
            value_name='Value'
        )

        melted_df_C4_only = pd.melt(
            selected_C4_only,
            id_vars=['V_BS_M', 'Photosynthesis'],
            value_vars=['Light Mesophyll', 'Light Bundle Sheath'],
            var_name='Light_Cell',
            value_name='Value'
        )

        melted_df_C3["Classification"] = "Only C3"
        melted_df_C4["Classification"] = "Variable"
        melted_df_C4_only["Classification"] = "Only C4"

        df = pd.concat([melted_df_C3, melted_df_C4, melted_df_C4_only])

        symbol_mapping = {
            'C3': 'circle',
            'Unidentified': 'cross',
            'C3_C4_Type_I': 'square',
            'C3_C4_Type_II': 'diamond',
            'C4': 'triangle-up'
        }

        # Create a new column in the DataFrame that maps Photosynthesis to symbols
        df['symbol'] = df['Photosynthesis'].map(symbol_mapping)


        # Creating an adjusted ATP value for the point chart
        df['Adjusted_Value'] = df['Value'] + (df['Value'].max() * 0.03)

        # Creating the bar chart
        bar = alt.Chart(df).mark_bar().encode(
            x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
            y=alt.Y('sum(Value):Q', axis=alt.Axis(grid=False, title="Photon Uptake (µmol m-2 s-1)"),
                   scale=alt.Scale(domain=[min_y, max_y])),
            color=alt.Color('Light_Cell:N', scale=alt.Scale(range=['#B8860B', '#2F4F4F']))
        )

        # Creating the point chart with a custom shape legend
        symbols = alt.Chart(df).mark_point(size=100).encode(
            x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
            y=alt.Y('sum(Adjusted_Value):Q', axis=alt.Axis(grid=False),
                   scale=alt.Scale(domain=[min_y, max_y])
                   ),
            shape=alt.Shape(
                'Photosynthesis:N',  # Use the Photosynthesis column to map to the shapes
                scale=alt.Scale(
                    domain=['C3', 'Unidentified', 'C3_C4_Type_I', 'C3_C4_Type_II', 'C4'],
                    range=['circle', 'cross', 'square', 'diamond', 'triangle-up']
                ),
                legend=alt.Legend(title="Photosynthesis Type")  # Custom legend title
            )
        )

        # Combining the bar and point charts using facet to ensure a shared y-axis
        chart = alt.layer(bar, symbols).facet(
            column=alt.Column('V_BS_M:N', title="Volume Ratio")
        ).resolve_scale(
            y='shared'  # Ensures that the y-axis is the same across all columns
        ).configure_view(
            strokeOpacity=0  # Remove grid lines around column clusters
        )

        return chart

def comparison_SoF_barplots_3D(df_C3, df_C4, df_C4_only ,PPFD, Vc_Vo, CO2, PC_BS, NH3, min_y, max_y):
        
        volume = [0.5, 1, 1.5, 2]
        
        selected_C3 = select_data(prepared_df=df_C3, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
        selected_C4 = select_data(prepared_df=df_C4, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
        selected_C4_only = select_data(prepared_df=df_C4_only, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)

        selected_C3 = selected_C3[selected_C3.V_BS_M.isin(volume)]
        selected_C4 = selected_C4[selected_C4.V_BS_M.isin(volume)]
        selected_C4_only = selected_C4_only[selected_C4_only.V_BS_M.isin(volume)]


        # Use melt to reshape the data
        melted_df_C4 = pd.melt(
            selected_C4,
            id_vars=['V_BS_M', 'Photosynthesis'],
            value_vars=[ "SoF M", "SoF BS","SoF transport"],
            var_name='SoF',
            value_name='Value'
        )

        melted_df_C3 = pd.melt(
            selected_C3,
            id_vars=['V_BS_M', 'Photosynthesis'],
            value_vars=[ "SoF M", "SoF BS","SoF transport"],
            var_name='SoF',
            value_name='Value'
        )

        melted_df_C4_only = pd.melt(
            selected_C4_only,
            id_vars=['V_BS_M', 'Photosynthesis'],
            value_vars=[ "SoF M", "SoF BS","SoF transport"],
            var_name='SoF',
            value_name='Value'
        )

        melted_df_C3["Classification"] = "Only C3"
        melted_df_C4["Classification"] = "Variable"
        melted_df_C4_only["Classification"] = "Only C4" 

        df = pd.concat([melted_df_C3, melted_df_C4, melted_df_C4_only])

        symbol_mapping = {
            'C3': 'circle',
            'Unidentified': 'cross',
            'C3_C4_Type_I': 'square',
            'C3_C4_Type_II': 'diamond',
            'C4': 'triangle-up'
        }

        # Create a new column in the DataFrame that maps Photosynthesis to symbols
        df['symbol'] = df['Photosynthesis'].map(symbol_mapping)


        # Creating an adjusted ATP value for the point chart
        df['Adjusted_Value'] = df['Value'] + (df['Value'].max() * 0.03)

        # Creating the bar chart
        bar = alt.Chart(df).mark_bar().encode(
            x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
            y=alt.Y('sum(Value):Q', axis=alt.Axis(grid=False, title="Sum of Fluxes (µmol m-2 s-1)"),
                   scale=alt.Scale(domain=[min_y, max_y])),
            color=alt.Color('SoF:N', scale=alt.Scale(range=['#D95F4D', '#5DADEC', "#F4A261"]))
        )

        # Creating the point chart with a custom shape legend
        symbols = alt.Chart(df).mark_point(size=100).encode(
            x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
            y=alt.Y('sum(Adjusted_Value):Q', axis=alt.Axis(grid=False),
                   scale=alt.Scale(domain=[min_y, max_y])),
            shape=alt.Shape(
                'Photosynthesis:N',  # Use the Photosynthesis column to map to the shapes
                scale=alt.Scale(
                    domain=['C3', 'Unidentified', 'C3_C4_Type_I', 'C3_C4_Type_II', 'C4'],
                    range=['circle', 'cross', 'square', 'diamond', 'triangle-up']
                ),
                legend=alt.Legend(title="Photosynthesis Type")  # Custom legend title
            )
        )

        # Combining the bar and point charts using facet to ensure a shared y-axis
        chart = alt.layer(bar, symbols).facet(
            column=alt.Column('V_BS_M:N', title="Volume Ratio")
        ).resolve_scale(
            y='shared'  # Ensures that the y-axis is the same across all columns
        ).configure_view(
            strokeOpacity=0  # Remove grid lines around column clusters
        )

        return chart

def reaction_lineplot_Vc_Vo(data, PPFD, V_BS_M, CO2, PC_BS, NH3):
    """ 
    Possible values for selected parameters (initial parameter scans):

    PPFD = [100, 200, 300, 400, 500, 600, 700, 800, 1000]
    CO2 = [20, 1000]
    Vc_Vo = [50, 55, 60 ,65, 70, 75, 80, 85, 90, 95, 99, 100]
    PC_BS = [1, .1]
    """
        
    selected_data = select_data_Vc_Vo(prepared_df=data, PPFD=PPFD, V_BS_M=V_BS_M, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)

    VM_list = []
    VBS_list = []

    for volume_ratio in selected_data.V_BS_M:
            VM_list.append(1/(1+volume_ratio))
            VBS_list.append(volume_ratio/(1+volume_ratio))

    #Create BS carboxylation and oxygenation variables
    selected_data['BS_Carboxylation'] = selected_data['BS_Rubisco'] + selected_data['BS_Rubisco_C4']
    selected_data['M_Carboxylation'] = selected_data['M_Rubisco'] 
    selected_data['BS_Oxygenation'] = selected_data['BS Oxygenase']
    selected_data['M_Oxygenation'] = selected_data['M Oxygenase']



    for n in range(len(selected_data.V_BS_M)):
            selected_data['M_Rubisco'].iloc[n] = selected_data['M_Rubisco'].iloc[n] * VM_list[n] 
            selected_data['BS_Rubisco'].iloc[n] = selected_data['BS_Rubisco'].iloc[n] * VBS_list[n] 
            selected_data['BS_Rubisco_C4'].iloc[n] = selected_data['BS_Rubisco_C4'].iloc[n] *  VBS_list[n]
            selected_data['Glycine Decarboxylase'].iloc[n] = selected_data['Glycine Decarboxylase'].iloc[n]  * VBS_list[n] 
            selected_data['PEP-CK'].iloc[n] = selected_data['PEP-CK'].iloc[n]  * VBS_list[n] 
            selected_data['NADP-ME'].iloc[n] = selected_data['NADP-ME'].iloc[n]  * VBS_list[n] 
            selected_data['NAD-ME'].iloc[n] = selected_data['NAD-ME'].iloc[n] * VBS_list[n]
            selected_data['BS_Carboxylation'].iloc[n] = selected_data['BS_Carboxylation'].iloc[n] * VBS_list[n]
            selected_data['M_Carboxylation'].iloc[n] = selected_data['M_Carboxylation'].iloc[n] * VM_list[n]
            selected_data['BS_Oxygenation'].iloc[n] = selected_data['BS_Oxygenation'].iloc[n] * VBS_list[n]
            selected_data['M_Oxygenation'].iloc[n] = selected_data['M_Oxygenation'].iloc[n] * VM_list[n]

    photo_df_melted = pd.melt(selected_data, id_vars=['Vc_Vo'], value_vars=[ 'M_Carboxylation', 'M_Oxygenation', 'BS_Carboxylation','BS_Oxygenation', "Carbon Uptake", "Glycine Decarboxylase", "PEP-CK", "NADP-ME", "NAD-ME"],
                    var_name='Rubisco_Type', value_name='Rubisco_Value')

    photo_df_melted.rename(columns={"Rubisco_Type": "Reactions"}, inplace=True)
    photo_df_melted["Reactions"] = photo_df_melted["Reactions"].map({"BS_Carboxylation": "BS Carboxylase", "M_Oxygenation": "M Oxygenase", "M_Carboxylation": "M Carboxylase", "BS_Oxygenation": "BS Oxygenase",  "Carbon Uptake" : "CO$_{2}$ Uptake", 
                                                                    "Glycine Decarboxylase": "Glycine Decarboxylase", "PEP-CK":"PEP-CK", "NADP-ME":"NADP-ME", "NAD-ME":"NAD-ME"})


    sns.set_style("ticks")

    line_styles = {"BS Carboxylase": "", "M Carboxylase": "", "BS Oxygenase":"", "M Oxygenase":"",  "CO$_{2}$ Uptake":"", 
                                                                "Glycine Decarboxylase": (2,2), "PEP-CK":(2,2), "NADP-ME":(2,2), "NAD-ME": (2,2)}
    
    # Define a custom color palette
    custom_palette = {"BS Carboxylase": "#bf5846", "M Carboxylase": "#30766a", "M Oxygenase": "#7eb2dd","BS Oxygenase": "#731727", "CO$_{2}$ Uptake": "#e42749",
                      "Glycine Decarboxylase": "#113630", "PEP-CK": "#5d679c", "NADP-ME": "#598322", "NAD-ME": "#5336c5"}

    g = sns.relplot(data=photo_df_melted, x="Vc_Vo", y="Rubisco_Value", hue="Reactions", kind="line", linewidth = 3, style="Reactions" ,dashes= line_styles, palette=custom_palette, alpha= .7)
    g.tight_layout()
    g.set_xticklabels(rotation = 45)
    g.set(xlabel="V$_{c}$/V$_{o}$ (%)", ylabel="Flux Value \n µmol cell type⁻¹ s$^{-1}$", xlim=(selected_data.Vc_Vo.min(), selected_data.Vc_Vo.max()), ylim=(0,  selected_data["Carbon Uptake"].max() + selected_data["Carbon Uptake"].max() * 1 ))

     
def reaction_lineplot_VRatio(data, PPFD, Vc_Vo, CO2, PC_BS, NH3):
    """ 
    Possible values for selected parameters (initial parameter scans):

    PPFD = [100, 200, 300, 400, 500, 600, 700, 800, 1000]
    CO2 = [20, 1000]
    Vc_Vo = [50, 55, 60 ,65, 70, 75, 80, 85, 90, 95, 99, 100]
    PC_BS = [1, .1]
    """
        
    selected_data = select_data(prepared_df=data, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3=NH3)
    
    VM_list = []
    VBS_list = []

    for volume_ratio in selected_data.V_BS_M:
            VM_list.append(1/(1+volume_ratio))
            VBS_list.append(volume_ratio/(1+volume_ratio))

    #Create BS carboxylation and oxygenation variables
    selected_data['BS_Carboxylation'] = selected_data['BS_Rubisco'] + selected_data['BS_Rubisco_C4']
    selected_data['M_Carboxylation'] = selected_data['M_Rubisco'] 
    selected_data['BS_Oxygenation'] = selected_data['BS Oxygenase']
    selected_data['M_Oxygenation'] = selected_data['M Oxygenase']

    for n in range(len(selected_data.V_BS_M)):
            selected_data['M_Rubisco'].iloc[n] = selected_data['M_Rubisco'].iloc[n] * VM_list[n] 
            selected_data['BS_Rubisco'].iloc[n] = selected_data['BS_Rubisco'].iloc[n] * VBS_list[n] 
            selected_data['BS_Rubisco_C4'].iloc[n] = selected_data['BS_Rubisco_C4'].iloc[n] *  VBS_list[n]
            selected_data['Glycine Decarboxylase'].iloc[n] = selected_data['Glycine Decarboxylase'].iloc[n]  * VBS_list[n] 
            selected_data['PEP-CK'].iloc[n] = selected_data['PEP-CK'].iloc[n]  * VBS_list[n] 
            selected_data['NADP-ME'].iloc[n] = selected_data['NADP-ME'].iloc[n]  * VBS_list[n] 
            selected_data['NAD-ME'].iloc[n] = selected_data['NAD-ME'].iloc[n] * VBS_list[n]
            selected_data['BS_Carboxylation'].iloc[n] = selected_data['BS_Carboxylation'].iloc[n] * VBS_list[n]
            selected_data['M_Carboxylation'].iloc[n] = selected_data['M_Carboxylation'].iloc[n] * VM_list[n]
            selected_data['BS_Oxygenation'].iloc[n] = selected_data['BS_Oxygenation'].iloc[n] * VBS_list[n]
            selected_data['M_Oxygenation'].iloc[n] = selected_data['M_Oxygenation'].iloc[n] * VM_list[n]
            selected_data['PEPC'].iloc[n] = selected_data['PEPC'].iloc[n] * VM_list[n]
    

    photo_df_melted = pd.melt(selected_data, id_vars=['V_BS_M'], value_vars=[ 'M_Carboxylation', 'M_Oxygenation', 'BS_Carboxylation','BS_Oxygenation', "Carbon Uptake", "Glycine Decarboxylase", "PEP-CK", "NADP-ME", "NAD-ME", "PEPC"],
                    var_name='Rubisco_Type', value_name='Rubisco_Value')

    photo_df_melted.rename(columns={"Rubisco_Type": "Reactions"}, inplace=True)
    photo_df_melted["Reactions"] = photo_df_melted["Reactions"].map({"BS_Carboxylation": "BS Carboxylase", "M_Oxygenation": "M Oxygenase", "M_Carboxylation": "M Carboxylase", "BS_Oxygenation": "BS Oxygenase",  "Carbon Uptake" : "CO$_{2}$ Uptake", 
                                                                    "Glycine Decarboxylase": "Glycine Decarboxylase", "PEP-CK":"PEP-CK", "NADP-ME":"NADP-ME", "NAD-ME":"NAD-ME", "PEPC":"PEPC"})


    sns.set_style("ticks")

    line_styles = {"BS Carboxylase": "", "M Carboxylase": "", "BS Oxygenase":"", "M Oxygenase":"",  "CO$_{2}$ Uptake":"", 
                                                                "Glycine Decarboxylase": (2,2), "PEP-CK":(2,2), "NADP-ME":(2,2), "NAD-ME": (2,2), "PEPC": (2,2)}
    
    # Define a custom color palette
    custom_palette = {"BS Carboxylase": "#bf5846", "M Carboxylase": "#30766a", "M Oxygenase": "#7eb2dd","BS Oxygenase": "#731727", "CO$_{2}$ Uptake": "#e42749",
                      "Glycine Decarboxylase": "#113630", "PEP-CK": "#5d679c", "NADP-ME": "#598322", "NAD-ME": "#5336c5", "PEPC": "#123456"}


    g = sns.relplot(data=photo_df_melted, x="V_BS_M", y="Rubisco_Value", hue="Reactions", kind="line", linewidth=3, style="Reactions", dashes=line_styles, palette=custom_palette, alpha= .7, height=4, aspect=5/4)    
    g.tight_layout()
    g.set_xticklabels(rotation = 45)
    g.set(xlabel="V$_{BS}$/V$_{M}$ (Unitless)", ylabel="Flux Value \n µmol cell type⁻¹ s$^{-1}$", xlim=(.1, 2), ylim=(0, 25))
    return g

def reaction_lineplot_VRatio_split(data, PPFD, Vc_Vo, CO2, PC_BS, NH3):
    selected_data = select_data(prepared_df=data, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3=NH3)
    selected_data = selected_data.copy()

    VM_list = []
    VBS_list = []

    for volume_ratio in selected_data.V_BS_M:
        VM_list.append(1/(1+volume_ratio))
        VBS_list.append(volume_ratio/(1+volume_ratio))

    # Create BS carboxylation and oxygenation variables
    selected_data['BS_Carboxylation'] = selected_data['BS_Rubisco'] + selected_data['BS_Rubisco_C4']
    selected_data['M_Carboxylation']  = selected_data['M_Rubisco']
    selected_data['BS_Oxygenation']   = selected_data['BS Oxygenase']
    selected_data['M_Oxygenation']    = selected_data['M Oxygenase']

    # Apply scaling safely using index labels
    for i, idx in enumerate(selected_data.index):
        selected_data.loc[idx, 'M_Rubisco']           *= VM_list[i]
        selected_data.loc[idx, 'BS_Rubisco']          *= VBS_list[i]
        selected_data.loc[idx, 'BS_Rubisco_C4']       *= VBS_list[i]
        selected_data.loc[idx, 'Glycine Decarboxylase'] *= VBS_list[i]
        selected_data.loc[idx, 'PEP-CK']               *= VBS_list[i]
        selected_data.loc[idx, 'NADP-ME']              *= VBS_list[i]
        selected_data.loc[idx, 'NAD-ME']               *= VBS_list[i]
        selected_data.loc[idx, 'BS_Carboxylation']     *= VBS_list[i]
        selected_data.loc[idx, 'M_Carboxylation']      *= VM_list[i]
        selected_data.loc[idx, 'BS_Oxygenation']       *= VBS_list[i]
        selected_data.loc[idx, 'M_Oxygenation']        *= VM_list[i]
        selected_data.loc[idx, 'PEPC']                 *= VM_list[i]

    # -------------------------------
    #       PLOT 1: RUBISCO
    # -------------------------------
    rubisco_vars = [
        'M_Carboxylation',
        'BS_Carboxylation',
        'BS_Oxygenation',
        'Carbon Uptake'
    ]

    df_rub = pd.melt(
        selected_data, id_vars=['V_BS_M'], value_vars=rubisco_vars,
        var_name='Reaction', value_name='Flux'
    )

    df_rub["Reaction"] = df_rub["Reaction"].map({
        "BS_Carboxylation": "BS Carboxylase",
        "M_Carboxylation": "M Carboxylase",
        "BS_Oxygenation": "BS Oxygenase",
        #"Carbon Uptake": "CO$_{2}$ Uptake"
    })

    # Strong contrasting palette
    rub_palette = {
        'BS Carboxylase': '#1f78b4',  # blue
        'M Carboxylase': '#33a02c',   # green
        'BS Oxygenase': '#ff7f00',    # orange
        #'CO$_{2}$ Uptake': 'red'
    }

    sns.set_style("ticks")

    g1 = sns.relplot(
        data=df_rub,
        x="V_BS_M", y="Flux", hue="Reaction",
        kind="line",
        linewidth=3,
        legend="brief",
        palette=rub_palette,
        markers=False,         # no markers
        dashes=False,           # solid lines
        height=4, 
        aspect=7.5/4)

    g1.tight_layout()
    g1.set(xticks=np.arange(0.2, 2.01, 0.2))
    g1.set_xticklabels(rotation=45)
    g1.set(
        xlabel="V$_{BS}$/V$_{M}$ (Unitless)",
        ylabel="Flux Value \n µmol cell type⁻¹ s$^{-1}$" ,
        xlim=(0.5, 2),
        ylim=(0, 25)
    )

    # -------------------------------
    #     PLOT 2: DECARBOXYLATION
    # -------------------------------
    decarb_vars = [
        "Glycine Decarboxylase",
        "PEP-CK",
        "NADP-ME",
        "NAD-ME",
        "PEPC",
        "Carbon Uptake"
    ]

    df_dec = pd.melt(
        selected_data, id_vars=['V_BS_M'], value_vars=decarb_vars,
        var_name='Reaction', value_name='Flux'
    )

    df_dec["Reaction"] = df_dec["Reaction"].map({
        "Glycine Decarboxylase": "Glycine Decarboxylase",
        "PEP-CK": "PEP-CK",
        "NADP-ME": "NADP-ME",
        "NAD-ME": "NAD-ME",
        "PEPC": "PEPC",
        #"Carbon Uptake": "CO$_{2}$ Uptake"
    })

    dec_palette = {
        'Glycine Decarboxylase': '#6a3d9a', # purple
        'PEP-CK': '#b15928',                # brown
        'NADP-ME': '#a6cee3',                # light blue
        'NAD-ME': '#fb9a99',                 # pink
        'PEPC': '#fdbf6f',                   # yellow-orange
        #'CO$_{2}$ Uptake': 'red'
    }

    g2 = sns.relplot(
        data=df_dec,
        x="V_BS_M", y="Flux", hue="Reaction",
        kind="line",
        linewidth=3,
        legend="brief",
        palette=dec_palette,
        markers=False,        
        dashes=False,           
        height=4, 
        aspect=7.5/4)          
    

    g2.tight_layout()
    g2.set(xticks=np.arange(0.2, 2.01, 0.2))
    g2.set_xticklabels(rotation=45)
    g2.set(
        xlabel="V$_{BS}$/V$_{M}$ (Unitless)",
        ylabel="Flux Value \n µmol cell type⁻¹ s$^{-1}$",
        xlim=(0.5, 2),
        ylim=(0, 25)
    )

    return g1, g2

def photosynthetic_bubble_plots(structured_C4, PPFD=1000, CO2=40, PC_BS=1, NH3=1000, 
                                desired_max_marker_px=18, output_svg="Bubble_Unsuberised.svg",
                                use_jitter=False):

    data = select_data_3D(prepared_df=structured_C4, PPFD=PPFD, CO2=CO2, PC_BS=PC_BS, NH3=NH3)
    
    # Filter out V_BS_M < 0.5
    data = data[data['V_BS_M'] >= 0.2].copy()

    VM_list = [1 / (1 + v) for v in data.V_BS_M]
    VBS_list = [v / (1 + v) for v in data.V_BS_M]
    data['BS_Carboxylation'] = data['BS_Rubisco'] + data['BS_Rubisco_C4']
    data['M_Carboxylation'] = data['M_Rubisco']
    data['BS_Oxygenation'] = data['BS Oxygenase']
    data['M_Oxygenation'] = data['M Oxygenase']

    for n in range(len(data)):
        data['M_Rubisco'].iloc[n] = data['M_Rubisco'].iloc[n] * VM_list[n]
        data['BS_Rubisco'].iloc[n] = data['BS_Rubisco'].iloc[n] * VBS_list[n]
        data['BS_Rubisco_C4'].iloc[n] = data['BS_Rubisco_C4'].iloc[n] * VBS_list[n]
        data['Glycine Decarboxylase'].iloc[n] = data['Glycine Decarboxylase'].iloc[n] * VBS_list[n]
        data['PEP-CK'].iloc[n] = data['PEP-CK'].iloc[n] * VBS_list[n]
        data['NADP-ME'].iloc[n] = data['NADP-ME'].iloc[n] * VBS_list[n]
        data['NAD-ME'].iloc[n] = data['NAD-ME'].iloc[n] * VBS_list[n]
        data['BS_Carboxylation'].iloc[n] = data['BS_Carboxylation'].iloc[n] * VBS_list[n]
        data['M_Carboxylation'].iloc[n] = data['M_Carboxylation'].iloc[n] * VM_list[n]
        data['BS_Oxygenation'].iloc[n] = data['BS_Oxygenation'].iloc[n] * VBS_list[n]
        data['M_Oxygenation'].iloc[n] = data['M_Oxygenation'].iloc[n] * VM_list[n]

    df = pd.DataFrame(data)
    photo_intermediates = ["C3", "C3_C4_Type_I", "C3_C4_Type_II", "Unidentified", "C4"]
    df_filtered = df.loc[df['Photosynthesis'].isin(photo_intermediates)]

    df_melted = df_filtered.melt(
        id_vars=['Vc_Vo', 'V_BS_M', "Photosynthesis"], 
        value_vars=['M_Carboxylation', 'BS_Carboxylation', "BS_Oxygenation"],
        var_name='Enzyme',
        value_name='Z'
    )

    df_melted = df_melted.rename(columns={"Z": "Flux Values"})
    df_melted["Enzyme"] = df_melted.Enzyme.replace({
        "M_Carboxylation": "Mesophyll Carboxylation", 
        "BS_Carboxylation": "Bundle Sheath Carboxylation", 
        "BS_Oxygenation": "Bundle Sheath Oxygenation"
    })
    df_melted["Photosynthesis"] = df_melted.Photosynthesis.replace({
        "C3_C4_Type_II": "C3-C4 Type II", 
        "C3_C4_Type_I": "C3-C4 Type I"
    })
    df_melted["Photosynthesis"] = df_melted.Photosynthesis.replace({
        "C3-C4 Type I": "C3-C4",
        "C3-C4 Type II": "C3-C4"
    })
    df_melted["Photosynthesis"] = pd.Categorical(
        df_melted["Photosynthesis"],
        categories=["C3", "Unidentified", "C3-C4", "C4"],
        ordered=True
    )

  
    if use_jitter:
        np.random.seed(0)
        jitter_mag = 0.01  #
        df_melted['V_BS_M_display'] = df_melted['V_BS_M'] + np.random.uniform(-jitter_mag, jitter_mag, size=len(df_melted))
        x_col = 'V_BS_M_display'
    else:
        x_col = 'V_BS_M'

    enzyme_colors = {
        'Mesophyll Carboxylation': '#06402B', 
        'Bundle Sheath Carboxylation': '#e6b705', 
        "Bundle Sheath Oxygenation": "#1f77b4"
    }
    custom_colors = ['#e6b705', "#06402B", '#1f77b4']

    # ---------- CONTROL MARKER SIZE (sizeref) ----------
    max_val = df_melted['Flux Values'].abs().max() if len(df_melted)>0 else 1.0
    desired_px = desired_max_marker_px
    sizeref = 3.0 * max_val / (desired_px ** 2 + 1e-12)

    # Create the faceted scatter
    fig = px.scatter(
        df_melted,
        y='Vc_Vo',
        x=x_col,
        size='Flux Values',
        size_max=desired_px,
        color='Enzyme',
        color_discrete_sequence=custom_colors,
        facet_col='Photosynthesis',
        category_orders={"Photosynthesis": ["C3", "Unidentified", "C3-C4", "C4"]},
        color_discrete_map=enzyme_colors,
        facet_col_spacing=0.01,   # increased spacing between facets
        height=420,
        width=1450  
    )

    # Apply sizeref and marker styling to all marker traces 
    fig.update_traces(
        marker=dict(sizemode='area', sizeref=sizeref, line=dict(width=0, color='black')),
        selector=dict(mode='markers')
    )

    # --- Axis formatting: extend x range more to prevent bubble clipping at edges ---
    fig.update_xaxes(range=[0.40, 2.1], dtick=0.1, tick0=0.5, tickformat=".1f", tickangle=-45)
    # ensure x title on every facet column 
    ncols = len(["C3", "Unidentified", "C3-C4", "C4"])
    for col in range(1, ncols + 1):
        fig.update_xaxes(title_text="<b>V<sub>BS/M</sub></b>", col=col)

    fig.update_layout(
        title=dict(text="", x=0.5, xanchor='center', font=dict(size=20, family='Arial')),
        yaxis_title="<b>V<sub>c</sub> / V<sub>o</sub></b>",
        xaxis_title="<b>V<sub>BS/M</sub></b>",
        legend_title=dict(text="<b>Enzymes</b>"),
        font=dict(size=14, family='Arial'),  
        plot_bgcolor='white',
        margin=dict(l=70, r=180, t=50, b=80),  # increased bottom margin for angled ticks
        legend=dict(
            orientation='v',  # 
            yanchor='middle',  #
            y=0.5,  
            x=1.02,  
            xanchor='left',  #
            font=dict(size=14, family='Arial')
        )
    )

    fig.for_each_annotation(lambda a: a.update(font=dict(size=16, family='Arial')))
    fig.update_xaxes(tickfont=dict(size=11, family='Arial'))
    fig.update_yaxes(tickfont=dict(size=11, family='Arial'))
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

 
    new_titles = {"C3": "<b>C<sub>3</sub></b>", "Unidentified": "<b>Unidentified</b>",
                  "C3-C4": "<b>C<sub>3</sub>-C<sub>4</sub></b>", "C4": "<b>C<sub>4</sub></b>"}
    fig.for_each_annotation(lambda a: a.update(text=new_titles.get(a.text, a.text)))

    fig.update_yaxes(mirror=True, ticks='outside', showline=True, linecolor='black',
                     gridcolor='#e8e8e8', griddash="dash", autorange="reversed",
                     title_font=dict(size=15, family='Arial'))
    fig.update_xaxes(mirror=True, ticks='outside', showline=True, linecolor='black',
                     gridcolor='#e8e8e8', griddash="dash",
                     title_font=dict(size=15, family='Arial'))

 
    expected_facets = ["C3", "Unidentified", "C3-C4", "C4"]
    present_facets = set(df_melted["Photosynthesis"].unique())
    for facet in expected_facets:
        if facet not in present_facets:
            col_num = expected_facets.index(facet) + 1
            fig.add_trace(
                go.Scatter(
                    y=[np.nan], x=[np.nan],
                    mode='markers',
                    marker=dict(opacity=0),
                    showlegend=False
                ), row=1, col=col_num
            )
    fig.update_xaxes(dtick=0.1)

    fig.update_xaxes(
    range=[0.45, 2.05],
    tickvals=np.arange(0.5, 2.01, 0.1),   
    ticktext=[f"{v:.1f}" for v in np.arange(0.5, 2.01, 0.1)],  
    tickangle=-45,
    tickformat=".1f"
    )
    

    # Add horizontal red line at Vc/Vo = 75 on all facets
    for col in range(1, len(expected_facets) + 1):
        fig.add_hline(
            y=75, 
            line_dash="solid", 
            line_color="red", 
            line_width=2,
            opacity=0.35, 
            row=1, 
            col=col,
            layer="below",
        )

    # Save as SVG -
    fig.write_image(output_svg, width=1500, height=420)
    fig.show()
    return fig



def comparison_ATP_split_barplots_3D_v2(df_C3, df_C4, PPFD, Vc_Vo, CO2, PC_BS, NH3, min_y, max_y):

    volume = [ 0.5, 1, 1.5, 2]
    
    selected_C3 = select_data(prepared_df=df_C3, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
    selected_C4 = select_data(prepared_df=df_C4, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)

    selected_C3 = selected_C3[selected_C3.V_BS_M.isin(volume)]
    selected_C4 = selected_C4[selected_C4.V_BS_M.isin(volume)]

    melted_df_C4 = pd.melt(
        selected_C4,
        id_vars=['V_BS_M', 'Photosynthesis'],
        value_vars=['ATP M', 'ATP BS'],
        var_name='ATP',
        value_name='Value'
    )

    melted_df_C3 = pd.melt(
        selected_C3,
        id_vars=['V_BS_M', 'Photosynthesis'],
        value_vars=['ATP M', 'ATP BS'],
        var_name='ATP',
        value_name='Value'
    )

    melted_df_C3["Classification"] = "Only C3"
    melted_df_C4["Classification"] = "Variable"

    df = pd.concat([melted_df_C3, melted_df_C4])

    # Merge C3-C4 types
    df['Photosynthesis'] = df['Photosynthesis'].replace({
        'C3_C4_Type_I': 'C3-C4',
        'C3_C4_Type_II': 'C3-C4'
    })

    bar = alt.Chart(df).mark_bar().encode(
        x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
        y=alt.Y('sum(Value):Q', axis=alt.Axis(grid=False, title="ATP Turnover (µmol m-2 s-1)"),
               scale=alt.Scale(domain=[min_y, max_y])
               ),
        color=alt.Color('Photosynthesis:N', 
                        scale=alt.Scale(
                            domain=['C3', 'Unidentified', 'C3-C4', 'C4'],
                            range=['#64a4ceff', '#884405fb', '#e6954cfb', '#69c569fe']
                        ),
                        legend=alt.Legend(title="Photosynthesis Type"))
    )

    chart = bar.facet(
        column=alt.Column('V_BS_M:N', title="Volume Ratio")
    ).resolve_scale(
        y='shared'  
    ).configure_view(
        strokeOpacity=0  
    )

    return chart

def comparison_Biomass_barplots_3D_v2(df_C3, df_C4, PPFD, Vc_Vo, CO2, PC_BS, NH3):
    """
    Possible values for selected parameters (initial parameter scans):

    PPFD = [100, 200, 300, 400, 500, 600, 700, 800, 1000]
    CO2 = [20, 1000]
    Vc_Vo = [50, 55, 60 ,65, 70, 75, 80, 85, 90, 95, 99, 100]
    PC_BS = [1, .1]

    """
    volume = [ 0.5, 1, 1.5, 2]
    
    selected_C3 = select_data(prepared_df=df_C3, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
    selected_C4 = select_data(prepared_df=df_C4, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)

    selected_C3 = selected_C3[selected_C3.V_BS_M.isin(volume)]
    selected_C4 = selected_C4[selected_C4.V_BS_M.isin(volume)]

    data_C4 = selected_C4.copy()
    data_C4.index=data_C4.V_BS_M
    data_C4["Classification"] = "Variable"

    data_C3 = selected_C3.copy()
    data_C3.index=data_C3.V_BS_M
    data_C3["Classification"] = "Only C3"

    df = pd.concat([data_C4, data_C3])

    # Merge C3-C4 types
    df['Photosynthesis'] = df['Photosynthesis'].replace({
        'C3_C4_Type_I': 'C3-C4',
        'C3_C4_Type_II': 'C3-C4'
    })

    # Creating the bar chart colored by Photosynthesis type
    bar = alt.Chart(df).mark_bar().encode(
        x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
        y=alt.Y('Biomass:Q', axis=alt.Axis(grid=False, title="Biomass Production (µmol m-2 s-1)")),
        color=alt.Color('Photosynthesis:N', 
                        scale=alt.Scale(
                            domain=['C3', 'Unidentified', 'C3-C4', 'C4'],
                            range=['#64a4ceff', '#884405fb', '#e6954cfb', '#69c569fe']
                        ),
                        legend=alt.Legend(title="Photosynthesis Type"))
    )

    # Combining using facet to ensure a shared y-axis
    chart = bar.facet(
        column=alt.Column('V_BS_M:N', title="Volume Ratio")
    ).resolve_scale(
        y='shared'  # Ensures that the y-axis is the same across all columns
    ).configure_view(
        strokeOpacity=0  # Remove grid lines around column clusters
    )

    return chart

def comparison_Light_barplots_3D_v2(df_C3, df_C4, PPFD, Vc_Vo, CO2, PC_BS, NH3, min_y, max_y):

        volume = [0.5, 1, 1.5, 2]
        
        selected_C3 = select_data(prepared_df=df_C3, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
        selected_C4 = select_data(prepared_df=df_C4, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)

        selected_C3 = selected_C3[selected_C3.V_BS_M.isin(volume)]
        selected_C4 = selected_C4[selected_C4.V_BS_M.isin(volume)]

        # Use melt to reshape the data
        melted_df_C4 = pd.melt(
            selected_C4,
            id_vars=['V_BS_M', 'Photosynthesis'],
            value_vars=['Light Mesophyll', 'Light Bundle Sheath'],
            var_name='Light_Cell',
            value_name='Value'
        )

        melted_df_C3 = pd.melt(
            selected_C3,
            id_vars=['V_BS_M', 'Photosynthesis'],
            value_vars=['Light Mesophyll', 'Light Bundle Sheath'],
            var_name='Light_Cell',
            value_name='Value'
        )

        melted_df_C3["Classification"] = "Only C3"
        melted_df_C4["Classification"] = "Variable"

        df = pd.concat([melted_df_C3, melted_df_C4])

        # Merge C3-C4 types
        df['Photosynthesis'] = df['Photosynthesis'].replace({
            'C3_C4_Type_I': 'C3-C4',
            'C3_C4_Type_II': 'C3-C4'
        })

        # Creating the bar chart colored by Photosynthesis type
        bar = alt.Chart(df).mark_bar().encode(
            x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
            y=alt.Y('sum(Value):Q', axis=alt.Axis(grid=False, title="Photon Uptake (µmol m-2 s-1)"),
                   scale=alt.Scale(domain=[min_y, max_y])),
            color=alt.Color('Photosynthesis:N', 
                            scale=alt.Scale(
                                domain=['C3', 'Unidentified', 'C3-C4', 'C4'],
                                range=['#64a4ceff', '#884405fb', '#e6954cfb', '#69c569fe']
                            ),
                            legend=alt.Legend(title="Photosynthesis Type"))
        )

        # Combining using facet to ensure a shared y-axis
        chart = bar.facet(
            column=alt.Column('V_BS_M:N', title="Volume Ratio")
        ).resolve_scale(
            y='shared'  # Ensures that the y-axis is the same across all columns
        ).configure_view(
            strokeOpacity=0  # Remove grid lines around column clusters
        )

        return chart

def comparison_SoF_barplots_3D_v2(df_C3, df_C4, PPFD, Vc_Vo, CO2, PC_BS, NH3, min_y, max_y):
        
        volume = [0.5, 1, 1.5, 2]
        
        selected_C3 = select_data(prepared_df=df_C3, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
        selected_C4 = select_data(prepared_df=df_C4, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)

        selected_C3 = selected_C3[selected_C3.V_BS_M.isin(volume)]
        selected_C4 = selected_C4[selected_C4.V_BS_M.isin(volume)]

        # Use melt to reshape the data
        melted_df_C4 = pd.melt(
            selected_C4,
            id_vars=['V_BS_M', 'Photosynthesis'],
            value_vars=[ "SoF M", "SoF BS","SoF transport"],
            var_name='SoF',
            value_name='Value'
        )

        melted_df_C3 = pd.melt(
            selected_C3,
            id_vars=['V_BS_M', 'Photosynthesis'],
            value_vars=[ "SoF M", "SoF BS","SoF transport"],
            var_name='SoF',
            value_name='Value'
        )

        melted_df_C3["Classification"] = "Only C3"
        melted_df_C4["Classification"] = "Variable"

        df = pd.concat([melted_df_C3, melted_df_C4])

        # Merge C3-C4 types
        df['Photosynthesis'] = df['Photosynthesis'].replace({
            'C3_C4_Type_I': 'C3-C4',
            'C3_C4_Type_II': 'C3-C4'
        })

        # Creating the bar chart colored by Photosynthesis type
        bar = alt.Chart(df).mark_bar().encode(
            x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
            y=alt.Y('sum(Value):Q', axis=alt.Axis(grid=False, title="Sum of Fluxes (µmol m-2 s-1)"),
                   scale=alt.Scale(domain=[min_y, max_y])),
            color=alt.Color('Photosynthesis:N', 
                            scale=alt.Scale(
                                domain=['C3', 'Unidentified', 'C3-C4', 'C4'],
                                range=['#64a4ceff', '#884405fb', '#e6954cfb', '#69c569fe']
                            ),
                            legend=alt.Legend(title="Photosynthesis Type"))
        )

        # Combining using facet to ensure a shared y-axis
        chart = bar.facet(
            column=alt.Column('V_BS_M:N', title="Volume Ratio")
        ).resolve_scale(
            y='shared'  # Ensures that the y-axis is the same across all columns
        ).configure_view(
            strokeOpacity=0  # Remove grid lines around column clusters
        )

        return chart

def comparison_QuantumYield_barplots_3D_v2(df_C3, df_C4, PPFD, Vc_Vo, CO2, PC_BS, NH3, min_y, max_y):

    volume = [0.5, 1, 1.5, 2]
    
    selected_C3 = select_data(prepared_df=df_C3, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)
    selected_C4 = select_data(prepared_df=df_C4, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3 = NH3)

    selected_C3 = selected_C3[selected_C3.V_BS_M.isin(volume)]
    selected_C4 = selected_C4[selected_C4.V_BS_M.isin(volume)]

    
    melted_df_C4 = pd.melt(
        selected_C4,
        id_vars=['V_BS_M', 'Photosynthesis'],
        value_vars=['Quantum Yield'],
        var_name='Metric',
        value_name='Value'
    )

    melted_df_C3 = pd.melt(
        selected_C3,
        id_vars=['V_BS_M', 'Photosynthesis'],
        value_vars=['Quantum Yield'],
        var_name='Metric',
        value_name='Value'
    )

    melted_df_C3["Classification"] = "Only C3"
    melted_df_C4["Classification"] = "Variable"

    df = pd.concat([melted_df_C3, melted_df_C4])

    # Merge C3-C4 types
    df['Photosynthesis'] = df['Photosynthesis'].replace({
        'C3_C4_Type_I': 'C3-C4',
        'C3_C4_Type_II': 'C3-C4'
    })

    # Creating the bar chart colored by Photosynthesis type
    bar = alt.Chart(df).mark_bar().encode(
        x=alt.X('Classification:N', title=None, axis=alt.Axis(labelAngle=-45)),
        y=alt.Y('Value:Q', axis=alt.Axis(grid=False, title="Quantum Yield"),
               scale=alt.Scale(domain=[min_y, max_y])
               ),
        color=alt.Color('Photosynthesis:N', 
                        scale=alt.Scale(
                            domain=['C3', 'Unidentified', 'C3-C4', 'C4'],
                            range=['#64a4ceff', '#884405fb', '#e6954cfb', '#69c569fe']
                        ),
                        legend=alt.Legend(title="Photosynthesis Type"))
    )

    # Combining using facet to ensure a shared y-axis
    chart = bar.facet(
        column=alt.Column('V_BS_M:N', title="Volume Ratio")
    ).resolve_scale(
        y='shared'  # Ensures that the y-axis is the same across all columns
    ).configure_view(
        strokeOpacity=0  # Remove grid lines around column clusters
    )

    return chart


def metabolic_comparison_subplots(df_C3, df_C4, PPFD, Vc_Vo, CO2, PC_BS, NH3):
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 20,
        'axes.labelsize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14})

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 8), sharex=False)
    volume_ratios = [0.5, 1, 1.5, 2]
    color_map = {'C3': '#64a4ceff', 'C3_C4_Type_I': '#e6954cfb', 'C4': '#69c569fe', 'Unidentified': '#884405fb'}
    
    def get_data():
        s_C3 = select_data(prepared_df=df_C3, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3=NH3)
        s_C4 = select_data(prepared_df=df_C4, PPFD=PPFD, Vc_Vo=Vc_Vo, CO2=CO2, PC_BS=PC_BS, NH3=NH3)
        s_C3 = s_C3[s_C3.V_BS_M.isin(volume_ratios)].copy()
        s_C4 = s_C4[s_C4.V_BS_M.isin(volume_ratios)].copy()
        s_C3['Classification'] = "C3"
        s_C4['Classification'] = "Variable"
        return pd.concat([s_C3, s_C4])

    df_combined = get_data()

    plot_configs = [
        (ax1, "Sum of Fluxes", "SoF transport", "SoF M", "SoF BS", True),
        (ax3, "Light Uptake", "Light Mesophyll", "Light Bundle Sheath", None, False),
        (ax2, "ATP Turnover", "ATP M", "ATP BS", None, False)]

    bar_width = 0.36  
    group_gap = 0.12  
    
    hatch_params = {'edgecolor': 'black', 'linewidth': 0}

    for ax, title, col1, col2, col3, has_tr in plot_configs:
        tick_positions = []
        tick_labels = []

        for i, vol in enumerate(volume_ratios):
            vol_data = df_combined[df_combined['V_BS_M'] == vol]
            for j, cls in enumerate(["C3", "Variable"]):
                row = vol_data[vol_data['Classification'] == cls]
                if row.empty: continue
                
                base_color = color_map.get(row['Photosynthesis'].iloc[0], 'black')
                pos = i * 1.2 + (j * (bar_width + group_gap)) 
                tick_positions.append(pos)
                tick_labels.append(cls)

                v1 = row[col1].iloc[0]
                v2 = row[col2].iloc[0]

                if has_tr: 
                    v3 = row[col3].iloc[0]
                    ax.bar(pos, v1, width=bar_width, color=base_color, hatch='..', **hatch_params)
                    ax.bar(pos, v2, bottom=v1, width=bar_width, color=base_color, hatch='//', **hatch_params)
                    ax.bar(pos, v3, bottom=v1+v2, width=bar_width, color=base_color, **hatch_params)
                    
                    # Draw the internal dividers
                    ax.hlines(y=[v1, v1 + v2], xmin=pos - bar_width/2, xmax=pos + bar_width/2, color='black', linewidth=1)
                else: 
                    ax.bar(pos, v1, width=bar_width, color=base_color, hatch='//', **hatch_params)
                    ax.bar(pos, v2, bottom=v1, width=bar_width, color=base_color, **hatch_params)
                    
                    # Draw the single internal divider
                    ax.hlines(y=v1, xmin=pos - bar_width/2, xmax=pos + bar_width/2, color='black', linewidth=1)

        ax.set_title(title, fontweight='bold', pad=25)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, rotation=45, ha='right')
        ax.spines[['top', 'right']].set_visible(False)

        for i, vol in enumerate(volume_ratios):
            group_center = i * 1.2 + (bar_width + group_gap) / 2
            ax.text(group_center, -0.18, f"{vol}", transform=ax.get_xaxis_transform(),
                    ha='center', va='top', fontsize=14)

    ax1.set_ylabel("Flux (µmol cell type⁻¹ s⁻¹)", fontsize=18)
    ax2.text(0.5, -0.3, r"$V_{BS/M}$ ratio", transform=ax2.transAxes, ha='center', fontsize=18)

    fill_elements = [
        Patch(facecolor='lightgray', edgecolor='black', linewidth=0.6, label='Bundle Sheath'),
        Patch(facecolor='lightgray', edgecolor='black', linewidth=0.6, hatch='//', label='Mesophyll'),
        Patch(facecolor='lightgray', edgecolor='black', linewidth=0.6, hatch='..', label='Transport')
    ]
    photo_elements = [Patch(facecolor=c, edgecolor='none', linewidth=0, label=l) for l, c in color_map.items()]
    
    ax3.legend(handles=photo_elements + [Patch(alpha=0, label="")] + fill_elements, 
                title="Photosynthetic Type", bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25, wspace=0.3, right=0.82)
    plt.savefig("../Images/metabolic_comp_subplots.svg", format="svg", bbox_inches="tight")
    plt.show()