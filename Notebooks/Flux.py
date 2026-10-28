import pandas as pd
import numpy as np
from typing import Optional


def calculate_ATP_turnover(df_solutions, metabolite, cell_type, c4_model):
    atp_turnover_list = []

    for column in df_solutions:
        # --- compute v_m, v_bs, pull out the per-reaction fluxes into df ---
        v_bs_m = float(column.split("_")[1])
        v_m    = 1.0 / (1.0 + v_bs_m)
        v_bs   = v_bs_m / (1.0 + v_bs_m)

        df = pd.DataFrame(df_solutions[column], index=df_solutions.index)

        # --- collect the metabolites of interest ---
        budget_metabolites = [
            met for met in c4_model.metabolites
            if metabolite in met.id and cell_type in met.id
        ]

        # --- build consumer/producer reaction lists (excluding transport) ---
        compartments = {"_c", "_m", "_p", "_x", "_e", "_i", "_v", "_l", "tx", "ss"}
        consumers, producers = [], []
        for met in budget_metabolites:
            for rxn in c4_model.reactions:
                tag = rxn.id[-2:]
                if tag in compartments:
                    if met in rxn.reactants:
                        consumers.append(rxn.id)
                    elif met in rxn.products:
                        producers.append(rxn.id)

        # --- pull out the fluxes ---
        prod_df = df.loc[producers].copy()
        cons_df = df.loc[consumers].copy()

        # any negative in prod → actually consumer, vice-versa
        neg_prods = prod_df[prod_df.iloc[:,0] < 0].index
        neg_cons  = cons_df[cons_df.iloc[:,0] < 0].index

        # flip them
        for rx in neg_prods:
            producers.remove(rx)
            consumers.append(rx)
        for rx in neg_cons:
            consumers.remove(rx)
            producers.append(rx)

        # re-slice, take abs, drop zero-flux rows
        prod_df = df.loc[producers].abs()
        prod_df = prod_df[(prod_df.T != 0).all()]

        cons_df = df.loc[consumers].abs()
        cons_df = cons_df[(cons_df.T != 0).all()]

        # tag them
        prod_df["Status"] = "Producer"
        cons_df["Status"] = "Consumer"

        # build master frame
        all_reactions = pd.concat([prod_df, cons_df], axis=0)
        if all_reactions.empty:
            # no reactions found at all → turnover is zero
            atp_turnover_list.append(0.0)
            continue

        # get stoichiometric coefficients
        labels = all_reactions.index
        atp_stoi = []
        for rx_id in labels:
            sto = 1.0
            try:
                rxn = c4_model.reactions.get_by_id(rx_id)
                # get_coefficient returns negative for reactant, positive for product
                sto = abs(rxn.get_coefficient(next(m for m in budget_metabolites if m in rxn.reactants + rxn.products)))
            except Exception:
                pass
            atp_stoi.append(sto)

        # compute the “new_flux” column
        base_flux = all_reactions.iloc[:, 0].values * atp_stoi
        #norm_factor = v_m if cell_type == "[M]" else v_bs
        all_reactions["new_flux"] = base_flux #* norm_factor

        # sum by Status, defaulting to 0 if missing
        grouped = all_reactions.groupby("Status")["new_flux"].sum()
        total_producers = grouped.get("Producer", 0.0)
        total_consumers = grouped.get("Consumer", 0.0)

        atp_turnover_list.append(total_producers)

    return atp_turnover_list

def get_sof(flux_solution_df):
    # Define lists of uptake and secretion reactions for mesophyll and bundle sheath
    M_Uptakes_List = ['[M]_Ca_tx', '[M]_H2O_tx', '[M]_CO2_tx', '[M]_O2_tx', '[M]_Pi_tx', 
                      '[M]_Mg_tx', '[M]_Nitrate_tx', '[M]_SO4_tx', '[M]_NH4_tx', '[M]_K_tx', 
                      '[M]_unlProtHYPO_c', '[M]_Photon_tx', '[M]_Arabidopsis_biomass_tx', '[M]_Maize_biomass_tx']

    BS_Uptakes_List = ['[B]_Ca_tx', '[B]_H2O_tx', '[B]_CO2_tx', '[B]_O2_tx', '[B]_Pi_tx', 
                       '[B]_Mg_tx', '[B]_Nitrate_tx', '[B]_SO4_tx', '[B]_NH4_tx', '[B]_K_tx', 
                       '[B]_unlProtHYPO_c', '[B]_Photon_tx']
    

    m_global_sof = []
    for column in flux_solution_df.columns[0:]:
        m_column_sof = 0
        for line in range(len(flux_solution_df.index)):
            if "[M]" in flux_solution_df.index[line] and flux_solution_df.index[line] not in M_Uptakes_List:
                m_column_sof += abs(flux_solution_df.loc[flux_solution_df.index[line]][column])
        m_global_sof.append(m_column_sof)

    bs_global_sof = []
    for column in flux_solution_df.columns[0:]:

        bs_column_sof = 0
        for line in range(len(flux_solution_df.index)):
            if "[B]" in flux_solution_df.index[line] and flux_solution_df.index[line] not in BS_Uptakes_List:
                bs_column_sof += abs(flux_solution_df.loc[flux_solution_df.index[line]][column])
        bs_global_sof.append(bs_column_sof)

    transport_global_sof = []

    for column in flux_solution_df.columns[0:]:
        t_column_sof = 0
        for line in range(len(flux_solution_df.index)):
            if "[MB]_" in flux_solution_df.index[line]:
                t_column_sof += abs(flux_solution_df.loc[flux_solution_df.index[line]][column])
        transport_global_sof.append(t_column_sof)

    return m_global_sof, bs_global_sof, transport_global_sof

def get_CET_LET(df_solutions, c4_model):

    cet_list = []
    let_list = []

    for column in df_solutions:

        try:

            df = pd.DataFrame(df_solutions[column], index=df_solutions.index)
            
            #Defining the list of metabolites to search
            budget_metabolites = []

            #Defining the compartments within the model
            compartment = ["_c", "_m", "_p", "_x", "_e", "_i", "_v", "_l", "tx", "ss"]

            #Searching for the metabolites in the model
            #for met in c4_model.metabolites.query("]_Oxidized_ferredoxins_p"):
            for met in c4_model.metabolites.query("]_PLASTOQUINOL_1_p"):
                if "[B]" in met.id:
                    budget_metabolites.append(met)


            #Defining list of reactions producing and consuming the metabolite
            consumers = []
            producers = []

            #Add reactions to respective list and exclude transport reactions
            for met in budget_metabolites:
                for reaction in c4_model.reactions:
                    if met in reaction.reactants and reaction.id[-2:] in compartment:
                        consumers.append(reaction.id)
                    elif met in reaction.products and reaction.id[-2:] in compartment:
                        producers.append(reaction.id)

            #Get flux values from the simulation for metabolite consuming/producing reactions
            producers_df =  df.loc[producers, :]
            #print(producers_df)
            consumers_df = df.loc[consumers, :]

            #Get values with negative flows: producing reactions with negative flow are consuming and vice-versa
            negative_producers = list(producers_df[producers_df.iloc[:, 0] < 0].index)
            negative_consumers = list(consumers_df[consumers_df.iloc[:, 0] < 0].index)

            #Add reactions to correct list
            consumers.extend(negative_producers)
            producers.extend(negative_consumers)

            #Remove reactions with negative flux from old list
            def remove_items(test_list, item):
                res = [i for i in test_list if i != item]
                return res

            for item in negative_producers:
                producers = remove_items(producers, item)

            for item in negative_consumers:
                consumers = remove_items(consumers, item)

            #Get flux values from the simulation for metabolite consuming/producing reactions (correct list)
            producers_df = df.loc[producers, :]
            #Make all values positive (disregard directionality)
            producers_df.iloc[:, 0] = producers_df.iloc[:, 0].abs()
            #Remove reactions with zero flux
            producers_df = producers_df[(producers_df.T != 0).all()]

            #Get flux values from the simulation for metabolite consuming/producing reactions (correct list)
            consumers_df = df.loc[consumers, :]
            #Make all values positive (disregard directionality)
            consumers_df.iloc[:, 0] = consumers_df.iloc[:, 0].abs()
            #Remove reactions with zero flux
            consumers_df = consumers_df[(consumers_df.T != 0).all()]

            producers_df["Status"] = "Producer"
            consumers_df["Status"] = "Consumer"

            frame = [producers_df, consumers_df]

            all_reactions = pd.concat(frame)

            all_reactions["label"] = all_reactions.index

            #Correct Budget stoichiometry
            atp_stoi = []
            for i in all_reactions['label']:
                for met in budget_metabolites:
                    try:
                        rxn = c4_model.reactions.get_by_id(i).get_coefficient(met)
                        atp_stoi.append(abs(rxn))
                    except KeyError:
                        continue

            new_flux = all_reactions.iloc[:, 0]* atp_stoi

            all_reactions.insert(3, 'coefficient', atp_stoi)
            all_reactions.insert(4, 'new_flux', new_flux)


            prod_reactions = all_reactions[all_reactions.Status == "Producer"]

            
            Total = all_reactions.groupby(["Status"]).new_flux.sum()[0]
            
            try:
                LET = prod_reactions._get_value('[B]_PSII_p', 'new_flux')
            except:
                LET = 0

            try:
                CET = prod_reactions._get_value('[B]_NADPH_Dehydrogenase_p', 'new_flux')
            except:
                CET = 0

            cet_list.append(CET)
            let_list.append(LET)
        
        except:
            cet_list.append(0)
            let_list.append(0)
    return cet_list, let_list

def get_PSI_PSII(df_solutions, c4_model):

    psi_list = []
    psii_list = []

    for column in df_solutions:

        try:

            df = pd.DataFrame(df_solutions[column], index=df_solutions.index)
            
            #Defining the list of metabolites to search
            budget_metabolites = []

            #Defining the compartments within the model
            compartment = ["_c", "_m", "_p", "_x", "_e", "_i", "_v", "_l", "tx", "ss"]

            #Searching for the metabolites in the model
            for met in c4_model.metabolites.query("]_Photon_p"):

                if "[B]" in met.id:
                    budget_metabolites.append(met)

            #Defining list of reactions producing and consuming the metabolite
            consumers = []
            producers = []

            #Add reactions to respective list and exclude transport reactions
            for met in budget_metabolites:
                for reaction in c4_model.reactions:
                    if met in reaction.reactants and reaction.id[-2:] in compartment:
                        consumers.append(reaction.id)
                    elif met in reaction.products and reaction.id[-2:] in compartment:
                        producers.append(reaction.id)

            #Get flux values from the simulation for metabolite consuming/producing reactions
            producers_df =  df.loc[producers, :]
            #print(producers_df)
            consumers_df = df.loc[consumers, :]

            #Get values with negative flows: producing reactions with negative flow are consuming and vice-versa
            negative_producers = list(producers_df[producers_df.iloc[:, 0] < 0].index)
            negative_consumers = list(consumers_df[consumers_df.iloc[:, 0] < 0].index)

            #Add reactions to correct list
            consumers.extend(negative_producers)
            producers.extend(negative_consumers)

            #Remove reactions with negative flux from old list
            def remove_items(test_list, item):
                res = [i for i in test_list if i != item]
                return res

            for item in negative_producers:
                producers = remove_items(producers, item)

            for item in negative_consumers:
                consumers = remove_items(consumers, item)

            #Get flux values from the simulation for metabolite consuming/producing reactions (correct list)
            producers_df = df.loc[producers, :]
            #Make all values positive (disregard directionality)
            producers_df.iloc[:, 0] = producers_df.iloc[:, 0].abs()
            #Remove reactions with zero flux
            producers_df = producers_df[(producers_df.T != 0).all()]

            #Get flux values from the simulation for metabolite consuming/producing reactions (correct list)
            consumers_df = df.loc[consumers, :]
            #Make all values positive (disregard directionality)
            consumers_df.iloc[:, 0] = consumers_df.iloc[:, 0].abs()
            #Remove reactions with zero flux
            consumers_df = consumers_df[(consumers_df.T != 0).all()]

            producers_df["Status"] = "Producer"
            consumers_df["Status"] = "Consumer"

            frame = [producers_df, consumers_df]

            all_reactions = pd.concat(frame)

            all_reactions["label"] = all_reactions.index

            #Correct Budget stoichiometry
            atp_stoi = []
            for i in all_reactions['label']:
                for met in budget_metabolites:
                    try:
                        rxn = c4_model.reactions.get_by_id(i).get_coefficient(met)
                        atp_stoi.append(abs(rxn))
                    except KeyError:
                        continue

            new_flux = all_reactions.iloc[:, 0]* atp_stoi

            all_reactions.insert(3, 'coefficient', atp_stoi)
            all_reactions.insert(4, 'new_flux', new_flux)


            consumer_reactions = all_reactions[all_reactions.Status == "Consumer"]

            
            Total = all_reactions.groupby(["Status"]).new_flux.sum()[0]

            try:
                PSII = consumer_reactions._get_value('[B]_PSII_RXN_p', 'new_flux')
            except:
                PSII = 0

            try:
                PSI = consumer_reactions._get_value('[B]_RXN490_3650_p', 'new_flux')
            except:
                PSI = 0

            psi_list.append(PSI)
            psii_list.append(PSII)
        
        except:
            psi_list.append(0)
            psii_list.append(0)
    
    return psi_list, psii_list

def ATP_Production(df, model):
    
    ATP_turn_M = calculate_ATP_turnover(cell_type = '[M]', df_solutions = df, metabolite = ']_ATP_', c4_model = model)
    aATP_turn_M = calculate_ATP_turnover(cell_type = '[M]', df_solutions = df, metabolite = ']_aATP_', c4_model = model)
    M_ATP_turnover = [x + y for x, y in zip(ATP_turn_M, aATP_turn_M)]

    ATP_turn_B = calculate_ATP_turnover(cell_type = '[B]', df_solutions = df, metabolite = ']_ATP_', c4_model = model)
    aATP_turn_B = calculate_ATP_turnover(cell_type = '[B]', df_solutions = df, metabolite = ']_aATP_', c4_model = model)
    B_ATP_turnover = [x + y for x, y in zip(ATP_turn_B, aATP_turn_B)]
    
    return M_ATP_turnover, B_ATP_turnover

def H_Production(df, model):
    
    NADPH_turn_M = calculate_ATP_turnover(cell_type = '[M]', df_solutions = df, metabolite = ']_NADPH_', c4_model = model)
    NADH_turn_M = calculate_ATP_turnover(cell_type = '[M]', df_solutions = df, metabolite = ']_NADH_', c4_model = model)
    M_H_turnover = [x + y for x, y in zip(NADPH_turn_M, NADH_turn_M)]

    NADPH_turn_BS = calculate_ATP_turnover(cell_type = '[B]', df_solutions = df, metabolite = ']_NADPH_', c4_model = model)
    NADH_turn_BS = calculate_ATP_turnover(cell_type = '[B]', df_solutions = df, metabolite = ']_NADH_', c4_model = model)
    BS_H_turnover = [x + y for x, y in zip(NADPH_turn_BS, NADH_turn_BS)]
    
    return M_H_turnover, BS_H_turnover

import pandas as pd


def photosynthetic_classification(df):
    f = df.copy()

    # Vectorized Filtering (GCV handled separately from [B])
    for col in f.columns:
        VM = 1 / (1 + float(col.split("_")[1]))
        VBS = float(col.split("_")[1]) / (1 + float(col.split("_")[1]))
        co2_10 = float(col.split("_")[5]) * 0.1
        co2_1 = float(col.split("_")[5]) * 0.05

        m_gcv = f.index == "[B]_GCVMULTI_RXN_m"
        m_b = f.index.str.contains(r"\[B\]") & ~m_gcv
        m_m = f.index.str.contains(r"\[M\]")

        f.loc[m_gcv & (f[col] * VBS).abs().le(co2_1), col] = 0
        f.loc[m_b & (f[col] * VBS).abs().le(co2_10), col] = 0
        f.loc[m_m & (f[col] * VM).abs().le(co2_10), col] = 0

    # Classification
    res = []
    for c in f.columns:
        # C3
        if (
            f.loc["[M]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p", c] > f.loc["[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_Ex", c]
            and f.loc["[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p", c] == 0
            and f.loc["[B]_PEPCARBOXYKIN_RXN_c", c] == 0
            and f.loc["[B]_MALIC_NADP_RXN_p", c] == 0
            and f.loc["[B]_1_PERIOD_1_PERIOD_1_PERIOD_39_RXN_m", c] == 0
            and f.loc["[B]_GCVMULTI_RXN_m", c] == 0
        ):
            res.append("C3")

        # C3-C4 Type I
        elif (
            f.loc["[M]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p", c] != 0
            and f.loc["[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p", c] != 0
            and f.loc["[B]_GCVMULTI_RXN_m", c] != 0
            and f.loc["[B]_PEPCARBOXYKIN_RXN_c", c] == 0
            and f.loc["[B]_MALIC_NADP_RXN_p", c] == 0
            and f.loc["[B]_1_PERIOD_1_PERIOD_1_PERIOD_39_RXN_m", c] == 0
        ):
            res.append("C3_C4_Type_I")

        # C3-C4 Type II
        elif (
            f.loc["[M]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p", c] != 0
            and f.loc["[M]_PEPCARBOX_RXN_c", c] != 0
            and f.loc["[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p", c] != 0
            and (
                f.loc["[B]_PEPCARBOXYKIN_RXN_c", c] != 0
                or f.loc["[B]_MALIC_NADP_RXN_p", c] != 0
                or f.loc["[B]_1_PERIOD_1_PERIOD_1_PERIOD_39_RXN_m", c] != 0
            )
            and f.loc["[B]_GCVMULTI_RXN_m", c] != 0
        ):
            res.append("C3_C4_Type_II")

        # C4
        elif (
            f.loc["[M]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p", c] == 0
            and f.loc["[M]_PEPCARBOX_RXN_c", c] != 0
            and f.loc["[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p", c] != 0
            and (
                f.loc["[B]_PEPCARBOXYKIN_RXN_c", c] != 0
                or f.loc["[B]_MALIC_NADP_RXN_p", c] != 0
                or f.loc["[B]_1_PERIOD_1_PERIOD_1_PERIOD_39_RXN_m", c] != 0
            )
        ):
            res.append("C4")

        # Unidentified
        else:
            res.append("Unidentified")

    return res
    
def photosynthetic_classification_old(df):   
    
    parameter_scan_filtered = df.copy()

    #Filtering as a percentage of CO2 assimilation (accounting for volume normalisation in the BS)
    for column in parameter_scan_filtered:
        CO2_filter = float(column.split("_")[5]) * 0.1 # filter by percentage of CO2 assimilation
        VM = 1/(1+float(column.split("_")[1]))
        VBS = float(column.split("_")[1])/(1+float(column.split("_")[1]))
        for line in parameter_scan_filtered[column].index:
            if "[B]" in line:
                if -abs(CO2_filter) <= parameter_scan_filtered[column].loc[line] *  VBS <= CO2_filter:
                    parameter_scan_filtered.at[line, column] = 0
                else:
                    continue
            elif "[M]" in line:
                if -abs(CO2_filter) <= parameter_scan_filtered[column].loc[line] *  VM <= CO2_filter:
                    parameter_scan_filtered.at[line, column] = 0

                else:
                    continue

    #Classify Flux Solutions
    #Preparing lists to create the dataset
    photo_classification = []

    for column in parameter_scan_filtered.columns[0:]:
        # Classifying the C3 - #Activity of Rubisco only in the Mesophyll 
        if (parameter_scan_filtered.loc["[M]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"][column] != 0 or parameter_scan_filtered.loc["[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_Ex"][column] != 0) and parameter_scan_filtered.loc["[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"][column] == 0 and parameter_scan_filtered.loc["[B]_PEPCARBOXYKIN_RXN_c"][column] == 0 and parameter_scan_filtered.loc["[B]_MALIC_NADP_RXN_p"][column] == 0 and parameter_scan_filtered.loc["[B]_1_PERIOD_1_PERIOD_1_PERIOD_39_RXN_m"][column] == 0 and parameter_scan_filtered.loc["[B]_GCVMULTI_RXN_m"][column] == 0:
            photo_classification.append("C3")
        # Classifying the C3-C4 type I - #Activity of M and BS Rubisco with operational Glycine Shuttle
        elif parameter_scan_filtered.loc["[M]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"][column] != 0 and parameter_scan_filtered.loc["[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"][column] != 0 and parameter_scan_filtered.loc["[B]_GCVMULTI_RXN_m"][column] != 0 and parameter_scan_filtered.loc["[B]_PEPCARBOXYKIN_RXN_c"][column] == 0 and parameter_scan_filtered.loc["[B]_MALIC_NADP_RXN_p"][column] == 0 and parameter_scan_filtered.loc["[B]_1_PERIOD_1_PERIOD_1_PERIOD_39_RXN_m"][column] == 0:
            photo_classification.append("C3_C4_Type_I")
        # Classifying the C3-C4 type II - #Activity of M and BS Rubisco with functional C4 cycle
        elif parameter_scan_filtered.loc["[M]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"][column] != 0 and parameter_scan_filtered.loc["[M]_PEPCARBOX_RXN_c"][column] != 0 and parameter_scan_filtered.loc["[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"][column] != 0 and (parameter_scan_filtered.loc["[B]_PEPCARBOXYKIN_RXN_c"][column] != 0 or parameter_scan_filtered.loc["[B]_MALIC_NADP_RXN_p"][column] != 0 or parameter_scan_filtered.loc["[B]_1_PERIOD_1_PERIOD_1_PERIOD_39_RXN_m"][column] != 0) and parameter_scan_filtered.loc["[B]_GCVMULTI_RXN_m"][column] != 0:
            photo_classification.append("C3_C4_Type_II")
        # Classifying the C4 - # Only BS Rubisco is operation - Exclusive activity of the C4 cycle
        elif parameter_scan_filtered.loc["[M]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"][column] == 0 and parameter_scan_filtered.loc["[M]_PEPCARBOX_RXN_c"][column] != 0 and parameter_scan_filtered.loc["[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"][column] != 0 and (parameter_scan_filtered.loc["[B]_PEPCARBOXYKIN_RXN_c"][column] != 0 or parameter_scan_filtered.loc["[B]_MALIC_NADP_RXN_p"][column] != 0 or parameter_scan_filtered.loc["[B]_1_PERIOD_1_PERIOD_1_PERIOD_39_RXN_m"][column] != 0):
            photo_classification.append("C4")
            # Unclassified flux solutions - Ex: Any situation in which CO2 is exclusively shuttled to the BS at the cost of photorespiration ("BS C3")
        else:
            photo_classification.append("Unidentified")
    
    return photo_classification

def structure_data(df: pd.DataFrame, model, filename: Optional[str] = None) -> pd.DataFrame:
    """
    Build the summary DataFrame from `df`. If `filename` is provided and points to a readable CSV,
    it is returned immediately; otherwise the summary is built and optionally saved to `filename`.
    """
    def _read_if_exists(fname):
        if fname is None:
            return None
        try:
            return pd.read_csv(fname, index_col=0)
        except Exception:
            return None

    def _prepare_indexed_copy(df_in: pd.DataFrame) -> pd.DataFrame:
        dfc = df_in.copy()
        if "Unnamed: 0" in dfc.columns:
            dfc = dfc.set_index("Unnamed: 0")
        return dfc

    def _parse_column_metadata(cols):
        # parse each column name into its components (expected at least 7 parts based on original code)
        # returns lists in the same order as columns
        parts = [col.split("_") for col in cols]
        # use safe indexing with fallback in case of malformed column names
        def safe_get(p, i, cast, default=np.nan):
            try:
                return cast(p[i])
            except Exception:
                return default
        light = [safe_get(p, 0, int) for p in parts]
        volume = [safe_get(p, 1, float) for p in parts]
        pc_m = [safe_get(p, 2, float) for p in parts]
        pc_bs = [safe_get(p, 3, float) for p in parts]
        photoresp = [safe_get(p, 4, float) for p in parts]
        max_c = [safe_get(p, 5, float) for p in parts]
        max_n = [safe_get(p, 6, float) for p in parts]
        return light, volume, pc_m, pc_bs, photoresp, max_c, max_n

    def _get_series_rounded(index_label, cols, decimals=1):
        # returns list of rounded values from df_copy.loc[index_label, col] for each col
        # falls back to np.nan if index label not present
        if index_label not in df_copy.index:
            return [np.nan] * len(cols)
        return [round(df_copy.at[index_label, col], decimals) for col in cols]

    def scale(x, in_min, in_max, out_min, out_max):
        return (x - in_min) / (in_max - in_min) * (out_max - out_min) + out_min

    #if filename exists and is readable, return it
    #existing = _read_if_exists(filename)
    #if existing is not None:
     #   return existing

    # build from df
    df_copy = _prepare_indexed_copy(df)
    cols = list(df_copy.columns)

    (light_max,
     volume_BS_M,
     PC_M_values,
     PC_BS_values,
     photoresp_ratio_values,
     maximum_carbon_assimilation,
     maximum_nitrogen_assimilation) = _parse_column_metadata(cols)

    # series extracted from df_copy with rounding
    light_mesophyll = _get_series_rounded("[M]_Photon_tx", cols, 1)
    light_bundle_sheath = _get_series_rounded("[B]_Photon_tx", cols, 1)
    carbon_uptakes = _get_series_rounded("[M]_CO2_tx", cols, 1)
    nitrogen_uptakes = _get_series_rounded("[B]_Nitrate_tx", cols, 1)
    biomass = _get_series_rounded("[M]_Arabidopsis_biomass_tx", cols, 4)

    Rubisco_M_list = _get_series_rounded("[M]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p", cols, 1)
    Rubisco_BS_C4_list = _get_series_rounded("[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p", cols, 1)
    Rubisco_BS_C3_list = _get_series_rounded("[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_Ex", cols, 1)
    PEP_CK_list = _get_series_rounded("[B]_PEPCARBOXYKIN_RXN_c", cols, 1)
    PEPC_list = _get_series_rounded("[M]_PEPCARBOX_RXN_c", cols, 1)
    NAD_ME_list = _get_series_rounded("[B]_1_PERIOD_1_PERIOD_1_PERIOD_39_RXN_m", cols, 1)
    NADP_ME_list = _get_series_rounded("[B]_MALIC_NADP_RXN_p", cols, 1)
    GLY_DECARBOX_list = _get_series_rounded("[B]_GCVMULTI_RXN_m", cols, 1)
    M_Oxygenase_list = _get_series_rounded("[M]_RXN_961_p", cols, 1)
    BS_Oxygenase_list = _get_series_rounded("[B]_RXN_961_p", cols, 1)

    # classification and other computed arrays
    photo_classification = photosynthetic_classification(df_copy)

    # scaling helpers and normalized volumes
    in_min, in_max = 0.1, 2
    out_min, out_max = 10, 1
    VM_Norm = [1 / (1 + v) for v in volume_BS_M]
    VBS_Norm = [v / (1 + v) for v in volume_BS_M]
    scaled_V_BS_M = [scale(v, in_min, in_max, out_min, out_max) for v in volume_BS_M]

    # other computed values using your domain functions
    m, bs, t = get_sof(df_copy)
    #cet, let = get_CET_LET(df_solutions=df_copy, c4_model=model)
    #psi, psii = get_PSI_PSII(df_solutions=df_copy, c4_model=model)

    cet = _get_series_rounded("[B]_NADPH_Dehydrogenase_p", cols, 1)
    let = _get_series_rounded("[B]_1_PERIOD_18_PERIOD_1_PERIOD_2_RXN_p", cols, 1)
    psi = _get_series_rounded("[B]_RXN490_3650_p", cols, 1)
    psii = _get_series_rounded("[B]_PSII_RXN_p", cols, 1)

    cet_2 = _get_series_rounded("[B]_NADPH_Dehydrogenase_p", cols, 1)
    let_2 = _get_series_rounded("[B]_1_PERIOD_18_PERIOD_1_PERIOD_2_RXN_p", cols, 1)


    ATP_turnover_M, ATP_turnover_BS = ATP_Production(df=df_copy, model=model)

    NADPH_turnover_M = _get_series_rounded("[M]_1_PERIOD_18_PERIOD_1_PERIOD_2_RXN_p", cols, 1)
    NADPH_turnover_BS = _get_series_rounded("[B]_1_PERIOD_18_PERIOD_1_PERIOD_2_RXN_p", cols, 1)

    Nitrate_reductase_B = _get_series_rounded("[B]_NITRATE_REDUCTASE_NADH_RXN_c", cols, 1)
    Nitrate_reductase_M = _get_series_rounded("[M]_NITRATE_REDUCTASE_NADH_RXN_c", cols, 1)
    Nitrite_reductase_B = _get_series_rounded("[B]_FERREDOXIN_NITRITE_REDUCTASE_RXN_p", cols, 1)
    Nitrite_reductase_M = _get_series_rounded("[M]_FERREDOXIN_NITRITE_REDUCTASE_RXN_p", cols, 1)

    Glutamine_syn_c_B = _get_series_rounded("[B]_GLUTAMINESYN_RXN_c", cols, 1)
    Glutamine_syn_m_B = _get_series_rounded("[B]_GLUTAMINESYN_RXN_m", cols, 1)
    Glutamine_syn_p_B = _get_series_rounded("[B]_GLUTAMINESYN_RXN_p", cols, 1)
    Glutamine_syn_c_M = _get_series_rounded("[M]_GLUTAMINESYN_RXN_c", cols, 1)
    Glutamine_syn_m_M = _get_series_rounded("[M]_GLUTAMINESYN_RXN_m", cols, 1)
    Glutamine_syn_p_M = _get_series_rounded("[M]_GLUTAMINESYN_RXN_p", cols, 1)

    glycolate_oxidase_M =  _get_series_rounded("[M]_RXN_969_x", cols, 1)
    glycolate_oxidase_B =  _get_series_rounded("[B]_RXN_969_x", cols, 1)
    phosphoglycolate_phosphatase_M = _get_series_rounded("[M]_GPH_RXN_p", cols, 1)
    phosphoglycolate_phosphatase_B = _get_series_rounded("[B]_GPH_RXN_p", cols, 1)

    # Build result DataFrame (only once)
    rows = list(zip(
        photoresp_ratio_values,
        light_max,
        volume_BS_M,
        maximum_nitrogen_assimilation,
        maximum_carbon_assimilation,
        carbon_uptakes,
        nitrogen_uptakes,
        np.array(light_mesophyll) * VM_Norm,
        np.array(light_bundle_sheath) * VBS_Norm,
        PC_M_values,
        PC_BS_values,
        np.array(m) * VM_Norm,
        np.array(bs) * VBS_Norm,
        np.array(t) * scaled_V_BS_M,
        photo_classification,
        biomass,
        np.array(ATP_turnover_BS) * VBS_Norm,
        np.array(ATP_turnover_M) * VM_Norm,
        np.array(NADPH_turnover_M) * VM_Norm, #
        np.array(NADPH_turnover_BS) * VBS_Norm, #
        Rubisco_M_list,
        Rubisco_BS_C3_list,
        Rubisco_BS_C4_list,
        GLY_DECARBOX_list,
        PEP_CK_list,
        NADP_ME_list,
        NAD_ME_list,
        PEPC_list,
        np.array(cet) * VBS_Norm,
        np.array(let) * VBS_Norm,
        np.array(psi) * VBS_Norm,
        np.array(psii) * VBS_Norm,
        BS_Oxygenase_list,
        M_Oxygenase_list,
        np.array(Nitrate_reductase_M) * VM_Norm, #
        np.array(Nitrate_reductase_B) * VBS_Norm, #
        np.array(Nitrite_reductase_M) * VM_Norm, #
        np.array(Nitrite_reductase_B) * VBS_Norm, #
        np.array(Glutamine_syn_c_M) * VM_Norm, #
        np.array(Glutamine_syn_m_M) * VM_Norm, #
        np.array(Glutamine_syn_p_M) * VM_Norm, #
        np.array(Glutamine_syn_c_B) * VBS_Norm, #
        np.array(Glutamine_syn_m_B) * VBS_Norm, #
        np.array(Glutamine_syn_p_B) * VBS_Norm, #
        glycolate_oxidase_M, 
        phosphoglycolate_phosphatase_M,
        glycolate_oxidase_B, 
        phosphoglycolate_phosphatase_B,
    ))

    columns = [
        "Vc_Vo", "PPFD", "V_BS_M", "Maximum Nitrogen Uptake", "Maximum Carbon Uptake",
        "Carbon Uptake", "Nitrogen Uptake", "Light Mesophyll", "Light Bundle Sheath",
        "PC_M", "PC_BS", "SoF M", "SoF BS", "SoF transport", "Photosynthesis", "Biomass",
        "ATP BS", "ATP M", "NADPH BS", "NADPH M", "M_Rubisco", "BS_Rubisco",
        "BS_Rubisco_C4", "Glycine Decarboxylase", "PEP-CK", "NADP-ME", "NAD-ME", "PEPC",
        "BS CET", "BS LET", "BS PSI", "BS PSII", "BS Oxygenase", "M Oxygenase",
        "Nitrate_reductase_M", "Nitrate_reductase_B", "Nitrite_reductase_M", "Nitrite_reductase_B", 
        "Glutamine_syn_c_M", "Glutamine_syn_m_M", "Glutamine_syn_p_M", "Glutamine_syn_c_B", "Glutamine_syn_m_B", "Glutamine_syn_p_B",
        "Glycolate oxidase M", "Phosphoglycolate phosphatase M", "Glycolate oxidase B", "Phosphoglycolate phosphatase B"
    ]

    result = pd.DataFrame(rows, columns=columns)

    # optionally save
    if filename is not None:
        result.to_csv(filename)

    return result


def select_data(prepared_df, PPFD, Vc_Vo, CO2, PC_BS, NH3): 
    if PPFD not in prepared_df["PPFD"].values:
        raise Exception("PPFD value is not present in df")
    
    if CO2 not in prepared_df["Maximum Carbon Uptake"].values: 
        raise Exception("CO2 uptake value is not present in df")
    
    if NH3 not in prepared_df["Maximum Nitrogen Uptake"].values: 
        raise Exception("Nitrogen Uptake value is not present in df")
    
    if Vc_Vo not in prepared_df["Vc_Vo"].values: 
        raise Exception("Photorespiratory ratio value is not present in df")
    
    if PC_BS not in prepared_df["PC_BS"].values: 
        raise Exception("Bundle sheath photosynthetic capacity value is not present in df")

    selected = prepared_df[
        (prepared_df["PPFD"] == PPFD) & 
        (prepared_df["Maximum Carbon Uptake"] == CO2) & 
        (prepared_df["Maximum Nitrogen Uptake"] == NH3) & 
        (prepared_df["Vc_Vo"] == Vc_Vo) & 
        (prepared_df["PC_BS"] == PC_BS)
    ]
    
    return selected


def select_data_Vc_Vo(prepared_df, PPFD, V_BS_M, CO2, PC_BS, NH3): 
    if PPFD not in prepared_df["PPFD"].values:
        raise Exception("PPFD value is not present in df")
    
    if NH3 not in prepared_df["Maximum Nitrogen Uptake"].values: 
        raise Exception("Nitrogen Uptake value is not present in df")
    
    if CO2 not in prepared_df["Maximum Carbon Uptake"].values: 
        raise Exception("CO2 uptake value is not present in df")
    
    if V_BS_M not in prepared_df["V_BS_M"].values:
        raise Exception("Volume ratio value is not present in df")
    
    if PC_BS not in prepared_df["PC_BS"].values:
        raise Exception("Bundle sheath photosynthetic capacity value is not present in df")

    selected = prepared_df[
        (prepared_df["PPFD"] == PPFD) & 
        (prepared_df["Maximum Carbon Uptake"] == CO2) &  
        (prepared_df["Maximum Nitrogen Uptake"] == NH3) & 
        (prepared_df["V_BS_M"] == V_BS_M) & 
        (prepared_df["PC_BS"] == PC_BS)
    ]

    return selected


def select_data_3D(prepared_df, PPFD, CO2, NH3, PC_BS): 
    if PPFD not in prepared_df["PPFD"].values:
        raise Exception("PPFD value is not present in df")
    
    if CO2 not in prepared_df["Maximum Carbon Uptake"].values: 
        raise Exception("CO2 uptake value is not present in df")

    if NH3 not in prepared_df["Maximum Nitrogen Uptake"].values: 
        raise Exception("Nitrogen Uptake value is not present in df")
    
    if PC_BS not in prepared_df["PC_BS"].values:
        raise Exception("Bundle sheath photosynthetic capacity value is not present in df")

    selected = prepared_df[
        (prepared_df["PPFD"] == PPFD) &  
        (prepared_df["Maximum Carbon Uptake"] == CO2) &  
        (prepared_df["Maximum Nitrogen Uptake"] == NH3) &  
        (prepared_df["PC_BS"] == PC_BS)
    ]

    return selected


def structure_data_transporters(df,filename): 
    
    #try:
     #   result = pd.read_csv(filename)
    
    #except:
    
    df_copy = df.copy()
    
    df_copy.index = df_copy["Unnamed: 0"]
    df_copy.drop(columns="Unnamed: 0", inplace = True)
    
    light_max = [int(column.split("_")[0]) for column in df_copy.columns[0:]]
    volume_BS_M = [float(column.split("_")[1]) for column in df_copy.columns[0:]]
    PC_M_values = [float(column.split("_")[2]) for column in df_copy.columns[0:]]
    PC_BS_values = [float(column.split("_")[3]) for column in df_copy.columns[0:]]
    photoresp_ratio_values = [float(column.split("_")[4]) for column in df_copy.columns[0:]]
    maximum_carbon_assimilation = [float(column.split("_")[5]) for column in df_copy.columns[0:]]
    maximum_nitrogen_assimilation = [float(column.split("_")[6]) for column in df_copy.columns[0:]]

    GLY_transport = [round(df_copy.loc["[MB]_GLY_c"][column], 1) for column in df_copy.columns[0:]]
    SER_transport = [round(df_copy.loc["[MB]_SER_c"][column], 1) for column in df_copy.columns[0:]]
    GAP_transport = [round(df_copy.loc["[MB]_GAP_c"][column], 1) for column in df_copy.columns[0:]]
    G3P_transport = [round(df_copy.loc["[MB]_G3P_c"][column], 1) for column in df_copy.columns[0:]]
    ALA_transport = [round(df_copy.loc["[MB]_L_ALPHA_ALANINE_c"][column], 1) for column in df_copy.columns[0:]]
    PYR_transport = [round(df_copy.loc["[MB]_PYRUVATE_c"][column], 1) for column in df_copy.columns[0:]]
    PEP_transport = [round(df_copy.loc["[MB]_PHOSPHO_ENOL_PYRUVATE_c"][column], 1) for column in df_copy.columns[0:]]
    ASP_transport = [round(df_copy.loc["[MB]_L_ASPARTATE_c"][column], 1) for column in df_copy.columns[0:]]
    MAL_transport = [round(df_copy.loc["[MB]_MAL_c"][column], 1) for column in df_copy.columns[0:]]
    GLYCOLLATE_transport = [round(df_copy.loc["[MB]_GLYCOLLATE_c"][column], 1) for column in df_copy.columns[0:]]


    result = pd.DataFrame(list(zip(photoresp_ratio_values, light_max, volume_BS_M, maximum_nitrogen_assimilation, maximum_carbon_assimilation, PC_M_values, PC_BS_values,
                                       GLY_transport, SER_transport, GAP_transport, G3P_transport, ALA_transport, PYR_transport, PEP_transport, ASP_transport, MAL_transport,GLYCOLLATE_transport)), 
                                        columns=["Vc_Vo", "PPFD", "V_BS_M", "Maximum Nitrogen Uptake" ,"Maximum Carbon Uptake", "PC_M", "PC_BS", "GLY", "SER", "GAP", "G3P", "ALA", "PYR", "PEP", "ASP", "MAL", "GLYCOL"])
    result.to_csv("{}".format(filename))

    return result


