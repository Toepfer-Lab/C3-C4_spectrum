"""
This module contains auxiliary functions to be used in Cobrapy model exploration
"""

def reaction_data(model, reaction_id, save = True, name = 'reaction_df.csv'):
    """
    This function receives a reaction ID and returns a Dataframe with the reaction metabolites information alongside
    writing it to a .csv file in the current working directory.

    :param
    model: A valid cobrapy model object
    reaction_id: An ID for a reaction in the model in string format (Ex : "Bio_opt")
    save: A boolean value to wheter the pandas dataframe should be saved to a .csv file. Default is True
    :return:
    A pandas Dataframe with the reaction's metabolites ID, Name, Compartment, Formula and Stoichiometric coefficient
    A "reaction_df.csv" file with the information in the pandas Dataframe
    """
    import pandas as pd
    import numpy as np
    import pathlib
    ID = []
    Name = []
    Compartment = []
    Formula = []
    Coefficient = []
    path = str(pathlib.Path().absolute()) + '\{}'.format(name)
    biomass = model.reactions.get_by_id(reaction_id)
    for met in biomass.metabolites:
        ID.append(met.id)
        Name.append(met.name)
        Compartment.append(met.compartment)
        Formula.append(met.formula)
        Coefficient.append(biomass.get_coefficient(met))
    df = pd.DataFrame(np.array([ID, Name, Compartment, Formula, Coefficient]))
    df = df.T
    df.columns = ['ID', 'Name', 'Compartment', 'Formula', 'Coefficient']
    if save:
        df.to_csv(r"{}".format(path), index = False, header = True)
    return df

def full_metabolite_data(model, save = True, name = 'metabolite_df.csv'):
    """
    This function receives a cobrapy model object and returns a Dataframe with the model metabolite information
    alongside writing it to a .csv file in the current working directory.

    :param
    model: A valid cobrapy model object
    save: A boolean value to wheter the pandas dataframe should be saved to a .csv file. Default is True.
    name: Name of the .csv file to store the dataframe. Default is 'metabolite_df.csv'
    :return:
    A pandas Dataframe with the model's metabolites ID, Name, Compartment and Formula
    A "metabolite_df.csv" file with the information in the pandas Dataframe
    """
    import pandas as pd
    import numpy as np
    import pathlib
    path = str(pathlib.Path().absolute()) + '\{}'.format(name)
    ID = []
    Name = []
    Compartment = []
    Formula = []
    for met in model.metabolites:
        ID.append(met.id)
        Name.append(met.name)
        Compartment.append(met.compartment)
        Formula.append(met.formula)
    df = pd.DataFrame(np.array([ID, Name, Compartment, Formula]))
    df = df.T
    df.columns = ['ID', 'Name', 'Compartment', 'Formula']
    if save:
        df.to_csv(r"{}".format(path), index = False, header = True)
    return df


def full_reaction_data(model, save = True, name = 'reactions_df.csv'):
    """

    @param model:
    @param save:
    @param name:
    @return:
    """
    import pandas as pd
    import numpy as np
    import pathlib
    path = str(pathlib.Path().absolute()) + '\{}'.format(name)
    ID = []
    Name = []
    Compartments = []
    Lower_Bound = []
    Upper_Bound = []
    for reaction in model.reactions:
        ID.append(reaction.id)
        Name.append(reaction.name)
        Compartments.append(reaction.compartments)
        Lower_Bound.append(reaction.lower_bound)
        Upper_Bound.append(reaction.upper_bound)
    df = pd.DataFrame(np.array([ID, Name, Compartments, Lower_Bound, Upper_Bound]))
    df = df.T
    df.columns = ['ID', 'Name',  'Compartments' , 'Lower_Bound', 'Upper_Bound']
    if save:
        df.to_csv(r"{}".format(path), index=False, header=True)
    return df

def test_from_metabolites(model, id1=None, id2=None):
    rec_list = []
    for reaction in model.reactions:
        for met in reaction.reactants:
            if id1 in met.id:
                rec_list.append(reaction)
    rec_list = list(dict.fromkeys(rec_list))
    if id2 != None:
        rec_list_2 = []
        for reaction in rec_list:
            for met in reaction.reactants:
                if id2 in met.id:
                    rec_list_2.append(reaction)
        rec_list_2 = list(dict.fromkeys(rec_list_2))
        for rec in rec_list_2:
            print(rec)
    else:
        for rec in rec_list:
            print(rec)

def set_fixed_flux(r_id, val, model):
    """
    This function receives a reaction ID, a numeric value and a Cobrapy model and fixes the bounds of the reaction
     corresponding to the reaction ID with the value provided
    @param r_id: Valid reaction id of the model being provided
    @param val: numeric value (Ex: 100)
    @param model: Valid Cobrapy model object
    """
    r_obj = model.reactions.get_by_id(r_id)
    r_obj.bounds = (val, val)


def set_bounds(r_id, val_tuple, model):
    """
    This function receives a reaction ID, a tuple with numeric values and a Cobrapy model and fixes the upper and lower
    bound of the reaction corresponding to the reaction ID with the values provided in the tuple
    @param r_id: Valid reaction id of the model being provided
    @param val_tuple: A tuple with numeric value (Ex: (-100,100))
    @param model:  Valid Cobrapy model object
    """
    r_obj = model.reactions.get_by_id(r_id)
    r_obj.bounds = val_tuple

def set_fixed_flux_ratio(r_dict, model):
    if len(r_dict) == 2:
        r_keys = r_dict.keys()
        r_values = r_dict.values()
        r_id1 = (list(r_keys))[0]
        r_obj1 = model.reactions.get_by_id(r_id1)
        r_v1 = (list(r_values))[0]
        r_id2 = (list(r_keys))[1]
        r_obj2 = model.reactions.get_by_id(r_id2)
        r_v2 = (list(r_values))[1]
        const = model.problem.Constraint(r_v1 * r_obj2.flux_expression - r_v2 * r_obj1.flux_expression, lb = 0, ub = 0)
        model.add_cons_vars(const)
        return const

def get_C4_fluxes(solution_frame):
    res = solution_frame.loc["[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"]["fluxes"]
    print(f"Bundle Sheath Rubisco carboxylase {res}")
    res = solution_frame.loc["[B]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_Ex"]["fluxes"]
    print(f"External Bundle Sheath Rubisco carboxylase {res}")
    res = solution_frame.loc["[M]_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"]["fluxes"]
    print(f"Mesophyll Rubisco carboxylase {res}")
    res = solution_frame.loc["[M]_PEPCARBOX_RXN_c"]["fluxes"]
    print(f"Mesophyll PEPC {res}")
    res = solution_frame.loc["[B]_PEPCARBOXYKIN_RXN_c"]["fluxes"]
    print(f"Bundle Sheath PEP-CK {res}")
    res = solution_frame.loc["[B]_MALIC_NADP_RXN_c"]["fluxes"]
    print(f"Bundle Sheath NADP-ME (c) {res}")
    res = solution_frame.loc["[B]_MALIC_NADP_RXN_p"]["fluxes"]
    print(f"Bundle Sheath NADP-ME (p) {res}")
    res = solution_frame.loc["[B]_1_PERIOD_1_PERIOD_1_PERIOD_39_RXN_m"]["fluxes"]
    print(f"Bundle Sheath NAD-ME {res}")
    res = solution_frame.loc["[B]_GCVMULTI_RXN_m"]["fluxes"]
    print(f"Bundle Sheath Glycine Dehydrogenase {res}")
    res = solution_frame.loc["[M]_GCVMULTI_RXN_m"]["fluxes"]
    print(f"Mesophyll Glycine Dehydrogenase {res}")


def check_maintenance_C4(solution_frame):
    res = solution_frame.loc["[B]_Photon_tx"]["fluxes"]
    print(f"Bundle Sheath photon uptake {res}")
    res = solution_frame.loc["[B]_ATPase_tx"]["fluxes"]
    print(f"Bundle Sheath ATP consumption {res}")
    res = solution_frame.loc["[B]_NADPHoxc_tx"]["fluxes"]
    print(f"Bundle Sheath cytoplasm NADPH consumption {res}")
    res = solution_frame.loc["[B]_NADPHoxp_tx"]["fluxes"]
    print(f"Bundle Sheath plastid NADPH consumption {res}")
    res = solution_frame.loc["[B]_NADPHoxm_tx"]["fluxes"]
    print(f"Bundle Sheath mitochondria NADPH consumption {res}")
    res = solution_frame.loc["[M]_Photon_tx"]["fluxes"]
    print(f"Mesophyll photon uptake {res}")
    res = solution_frame.loc["[M]_ATPase_tx"]["fluxes"]
    print(f"Mesophyll ATP consumption {res}")
    res = solution_frame.loc["[M]_NADPHoxc_tx"]["fluxes"]
    print(f"Mesophyll cytoplasm NADPH consumption {res}")
    res = solution_frame.loc["[M]_NADPHoxp_tx"]["fluxes"]
    print(f"Mesophyll plastid NADPH consumption {res}")
    res = solution_frame.loc["[M]_NADPHoxm_tx"]["fluxes"]
    print(f"Mesophyll mitochondria NADPH consumption {res}")

def check_maintenance_C3(solution_frame):
    res = solution_frame.loc["Photon_tx"]["fluxes"]
    print(f"Photon uptake {res}")
    res = solution_frame.loc["ATPase_tx"]["fluxes"]
    print(f"ATP consumption {res}")
    res = solution_frame.loc["NADPHoxc_tx"]["fluxes"]
    print(f"Cytoplasm NADPH consumption {res}")
    res = solution_frame.loc["NADPHoxp_tx"]["fluxes"]
    print(f"Plastid NADPH consumption {res}")
    res = solution_frame.loc["NADPHoxm_tx"]["fluxes"]
    print(f"Mitochondria NADPH consumption {res}")
