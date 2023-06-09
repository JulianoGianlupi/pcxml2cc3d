import warnings

# defines conversion factors to meter
_space_convs = {"micron": 1e-6,
                "micrometer": 1e-6,
                "micro": 1e-6,
                "milli": 1e-3,
                "millimeter": 1e-3,
                "nano": 1e-9,
                "nanometer": 1e-9,
                'meter': 1
                }

# defines conversion factors to minutes
_time_convs = {"millisecond": 1e-3 / 60,
               "milliseconds": 1e-3 / 60,
               "microsecond": 1e-6 / 60,
               "microseconds": 1e-6 / 60,
               "second": 1 / 60,
               "s": 1 / 60,
               "seconds": 1 / 60,
               "hours": 60,
               "hour": 60,
               "h": 60,
               "day": 24 * 60,
               "days": 24 * 60,
               "week": 7 * 24 * 60,
               "weeks": 7 * 24 * 60,
               "minutes": 1,
               "min": 1}


def get_dims(pcdict, space_convs=_space_convs):
    """
    Parses PhysiCell data and generates CC3D space dimensions and unit conversions

    This function looks for the value of the maximum and minimum of all coordinates in PhysiCell (
    pcdict['domain']['x_min'], pcdict['domain']['x_max'], etc) and saves them to variables. It also looks for the
    discretization variables from PhysiCell (pcdict['domain']['dx'], etc) and saves them to variables. Using the size of
    the domain and the discretization it defines what will be the number of pixels in CompuCell3D's domain.
    It also looks for the unit used in PhysiCell to determine what will be the pixel/unit factor in CC3D.

    :param pcdict: Dictionary created from parsing PhysiCell XML
    :param space_convs: Dictionary of predefined space units
    :return pcdims, ccdims: Two tuples representing the dimension data from PhysiCell and in CC3D.
          ((xmin, xmax), (ymin, ymax), (zmin, zmax), units), and
          (cc3dx, cc3dy, cc3dz, cc3dspaceunitstr, cc3dds, autoconvert_space, is_2D)
    """
    is_2D = True if pcdict['domain']['use_2D'].upper() == 'TRUE' else False
    xmin = float(pcdict['domain']['x_min']) if "x_min" in pcdict['domain'].keys() else None
    xmax = float(pcdict['domain']['x_max']) if "x_max" in pcdict['domain'].keys() else None

    ymin = float(pcdict['domain']['y_min']) if "y_min" in pcdict['domain'].keys() else None
    ymax = float(pcdict['domain']['y_max']) if "y_max" in pcdict['domain'].keys() else None

    zmin = float(pcdict['domain']['z_min']) if "z_min" in pcdict['domain'].keys() else None
    zmax = float(pcdict['domain']['z_max']) if "z_max" in pcdict['domain'].keys() else None

    units = pcdict['overall']['space_units'] if 'overall' in pcdict.keys() and \
                                                'space_units' in pcdict['overall'].keys() else 'micron'

    autoconvert_space = True
    if units not in space_convs.keys():
        message = f"WARNING: {units} is not part of known space units. Automatic space-unit conversion disabled." \
                  f"Available units for auto-conversion are:\n{space_convs.keys()}"
        warnings.warn(message)
        autoconvert_space = False

    # the dx/dy/dz tags mean that for every voxel there are dx space-units.
    # therefore [dx] = [space-unit/voxel]. Source: John Metzcar

    dx = float(pcdict['domain']['dx']) if "dx" in pcdict['domain'].keys() else 1
    dy = float(pcdict['domain']['dy']) if "dy" in pcdict['domain'].keys() else 1
    dz = float(pcdict['domain']['dz']) if "dz" in pcdict['domain'].keys() else 1

    # print(dx, dy, dz, type(dx), type(dy), type(dz), )
    if not dx == dy == dz:
        message = "WARNING! Physicell's dx/dy/dz are not all the same: " \
                  f"dx={dx}, dy={dy}, dz={dz}\n" \
                  f"Using {max(min(dx, dy, dz), 1)}"
        warnings.warn(message)

    ds = max(min(dx, dy, dz), 1)

    diffx = 1 if xmin is None or xmax is None else round(xmax - xmin)
    diffy = 1 if ymin is None or ymax is None else round(ymax - ymin)
    diffz = 1 if zmin is None or zmax is None else round(ymax - zmin)

    cc3dx = round(diffx / ds)
    cc3dy = round(diffy / ds)
    cc3dz = round(diffz / ds) if not is_2D else 1

    cc3dds = cc3dx / diffx  # pixel/unit

    # [cc3dds] = pixel/unit
    # [cc3dds] * unit = pixel

    cc3dspaceunitstr = f"1 pixel = {cc3dds} {units}"

    pcdims, ccdims = ((xmin, xmax), (ymin, ymax), (zmin, zmax), units), \
                     (cc3dx, cc3dy, cc3dz, cc3dspaceunitstr, cc3dds, autoconvert_space, is_2D)

    return pcdims, ccdims


def get_time(pcdict, time_convs=_time_convs):
    """
    Parses PhysiCell data and generates CC3D time dimensions and unit conversions

    This function fetches the maximum time set in PhysiCell (pcdict['overall']['max_time']['#text']) and the time
    discretization (pcdict['overall']['dt_mechanics']['#text']) to set the max time for the CC3D simulation and what is
    MCS/unit factor in CC3D

    :param pcdict: Dictionary created from parsing PhysiCell XML
    :param time_convs: Dictionary of predefined space units
    :return pctime, cctime: Two tuples representing the dimension data from PhysiCell and in CC3D.
        (mt, time_unit, mechdt), and (steps, cc3dtimeunitstr, cc3ddt, autoconvert_time)
    """

    mt = float(pcdict['overall']['max_time']['#text']) if "max_time" in pcdict['overall'].keys() and \
                                                          '#text' in pcdict['overall']['max_time'].keys() else 100000

    mtunit = pcdict['overall']['max_time']['@units'] if "max_time" in pcdict['overall'].keys() and '@units' in \
                                                        pcdict['overall']['max_time'].keys() else None

    time_unit = pcdict['overall']['time_units'] if "time_units" in pcdict['overall'].keys() else None

    if mtunit != time_unit:
        message = f"Warning: Psysicell time units in " \
                  "\n`<overall>\n\t<max_time units=...`\ndiffers from\n" \
                  f"`<time_units>unit</time_units>`.\nUsing: {time_unit}"
        warnings.warn(message)

    autoconvert_time = True
    if time_unit not in time_convs.keys():
        message = f"WARNING: {time_unit} is not part of known time units. Automatic time-unit conversion disabled." \
                  f"Available units for auto-conversion are:\n{time_convs.keys()}"
        warnings.warn(message)
        autoconvert_time = False

    mechdt = float(pcdict['overall']['dt_mechanics']['#text']) if "dt_mechanics" in pcdict['overall'].keys() else 0.1

    if autoconvert_time and time_unit != "min":
        mechdt *= time_convs[time_unit]

    steps = round(mt / mechdt)

    cc3ddt = 1 / (mt / steps)  # MCS/unit

    # [cc3ddt] = MCS/unit
    # [cc3ddt] * unit = MCS

    cc3dtimeunitstr = f"1 MCS = {cc3ddt} {time_unit}"

    # timeconvfact = 1/cc3ddt
    pctime, cctime = (mt, time_unit, mechdt), (steps, cc3dtimeunitstr, cc3ddt, autoconvert_time)

    return pctime, cctime


def get_parallel(pcdict):
    return int(pcdict['parallel']['omp_num_threads']) if 'parallel' in pcdict.keys() and \
                                                         'omp_num_threads' in pcdict['parallel'].keys() else 1


def get_boundary_wall(pcdict):
    """
    Looks for the existance of a wall around the PhysiCell simulation and returns a flag saying if it does or not
    :param pcdict: Dictionary created from parsing PhysiCell XML
    :return wall_exists: bool for the existance of the boundary wall
    """
    wall_exists = False
    if "options" in pcdict.keys() and "virtual_wall_at_domain_edge" in pcdict['options'].keys():
        wall = pcdict["options"]["virtual_wall_at_domain_edge"].upper()
        if wall == "TRUE":
            wall_exists = True
            return wall_exists
        else:
            return wall_exists
    else:
        return wall_exists
    return wall_exists


def get_cell_volume(subdict):
    if 'phenotype' in subdict.keys() and 'volume' in subdict['phenotype'].keys() and \
            'total' in subdict['phenotype']['volume'].keys():
        volume = float(subdict['phenotype']['volume']['total']['#text'])
        units = subdict['phenotype']['volume']['total']['@units']
        return volume, units
    return None, None


def get_cell_mechanics(subdict):
    """
    Extracts the mechanics data for a given cell from a pcdict['cell_definitions']['cell_definition'] subdictionary.

    Parameters:
    -----------
    subdict : dict
        A dictionary containing information about the cell.

    Returns:
    --------
    mechanics : dict or None
        A dictionary containing the mechanics data for the given cell, with keys representing the different mechanical
        properties (e.g. "cell_cell_adhesion_strength") and values representing the corresponding
        numerical values and units. Returns None if the given subdictionary does not contain mechanics data.
    """
    if 'phenotype' in subdict.keys() and 'mechanics' in subdict['phenotype'].keys():
        d = {}
        for key, item in subdict['phenotype']['mechanics'].items():
            if key != "options":
                d[key] = {'units': item['@units'],
                          'value': float(item['#text'])}
    else:
        return None

    return d


def check_below_minimum_volume(volume, minimum=8):
    """
        Checks if the given volume falls below the given minimum volume threshold.

        Parameters:
        -----------
        volume : float or None
            The volume to check. If None, returns False and the given minimum.
        minimum : float, optional
            The minimum volume threshold. If the given volume falls below this threshold, returns True and this value.
            Defaults to 8.

        Returns:
        --------
        below : bool
            A boolean indicating whether the given volume falls below the minimum volume threshold.
        minimum : float
            Returns the minimum volume.
    """
    if volume is None:
        return False, minimum
    return volume < minimum, minimum


_physicell_phenotype_codes = {
    "0": "Ki67 Advanced",
    "1": "Ki67 Basic",
    "2": "Flow Cytometry Basic",
    "5": "Simple Live",
    "6": "Flow Cytometry Advanced",
    "7": "Ki67 Basic",
    "100": "Standard apoptosis model",
    "101": "Standard necrosis model"}


def get_cycle_rate_data(rate_data, volume_datum, using_rates):
    """
    Returns a tuple of 10 lists containing data related to cycle rates, given input data.

    The get_cycle_rate_data function takes in three arguments: rate_data, volume_datum, and using_rates. The function
    returns a tuple of ten lists, which contain data related to cycle rates.

    The rate_data argument can be either a list of dictionaries or a dictionary. If it is a dictionary, the function
    calculates the phase duration using the value of the '#text' key and stores it in a variable phase_duration. If
    using_rates is True, phase_duration is calculated as the reciprocal of the '#text' value. Otherwise, phase_duration
    is calculated as the float value of the '#text' key. If '#text' is not present or if its value is 0, 9e99 is stored
    as the phase_duration.

    The fixed_duration key in the rate_data dictionary is used to store the fixed duration of the phase. The tuple
    (fixed_duration, phase_duration) is then stored in a list called phase_durations. Other relevant data from the
    volume_datum dictionary is also stored in separate lists.

    If rate_data is a list of dictionaries, the function iterates over each dictionary in the list and performs the
    same operations as above for each dictionary. The resulting data is then stored in separate lists.

    The function returns a tuple of ten lists: phase_durations, fluid_change_rate, cytoplasmic_biomass_change_rate,
    nuclear_biomass_change_rate, calcification_rate, fluid_fraction, nuclear, calcified_fraction, rel_rupture, and
    total. The docstring lists the data type and content of each list.

    Parameters:
    -----------
        rate_data : list/dict
            a list of dictionaries, or a dictionary, containing rate data.
        volume_datum : dict
            a dictionary containing volume data.
        using_rates : bool
            a boolean indicating whether rate data is being used.

    Returns:
    --------
    A tuple of 10 lists:

        phase_durations : list
            a list of tuples containing phase duration data
        fluid_change_rate : list
            a list of floats containing fluid change rate data.
        cytoplasmic_biomass_change_rate : list
            a list of floats containing cytoplasmic biomass change rate data.
        nuclear_biomass_change_rate : list
            a list of floats containing nuclear biomass change rate data.
        calcification_rate : list
            a list of floats containing calcification rate data.
        fluid_fraction : list
            a list of floats containing fluid fraction data.
        nuclear : list
            a list of floats containing nuclear data.
        calcified_fraction : list
            a list of floats containing calcified fraction data.
        rel_rupture : list
            a list of None values.
        total : list
            a list of floats containing total volume data.
    """
    if type(rate_data) != list:

        if using_rates:
            phase_duration = 1 / float(rate_data['#text']) if float(rate_data['#text']) \
                else 9e99
        else:
            phase_duration = float(rate_data['#text']) if float(rate_data['#text']) \
                else 9e99

        fixed_duration = rate_data['@fixed_duration'].upper()
        duration_data = (fixed_duration, phase_duration)
        phase_durations = [duration_data]
        fluid_change_rate = [float(volume_datum['fluid_change_rate']['#text'])]
        cytoplasmic_biomass_change_rate = [float(volume_datum['cytoplasmic_biomass_change_rate']['#text'])]
        nuclear_biomass_change_rate = [float(volume_datum['nuclear_biomass_change_rate']['#text'])]
        calcification_rate = [float(volume_datum['calcification_rate']['#text'])]
        fluid_fraction = [float(volume_datum['fluid_fraction']['#text'])]
        nuclear = [float(volume_datum['nuclear']['#text'])]
        calcified_fraction = [float(volume_datum['calcified_fraction']['#text'])]
        rel_rupture = [None]
        total = [float(volume_datum['total']['#text'])]
    else:

        phase_durations = []
        fluid_change_rate = []
        cytoplasmic_biomass_change_rate = []
        nuclear_biomass_change_rate = []
        calcification_rate = []

        fluid_fraction = []
        nuclear = []
        calcified_fraction = []
        rel_rupture = []
        total = []
        for rate_datum in rate_data:
            fixed_duration = rate_datum['@fixed_duration'].upper()
            if using_rates:
                phase_duration = 1 / float(rate_datum['#text']) if float(rate_datum['#text']) \
                    else 9e99
            else:
                phase_duration = float(rate_datum['#text']) if float(rate_datum['#text']) \
                    else 9e99
            duration_data = (fixed_duration, phase_duration)
            phase_durations.append(duration_data)

            fluid_change_rate.append(float(volume_datum['fluid_change_rate']['#text']))
            cytoplasmic_biomass_change_rate.append(
                float(volume_datum['cytoplasmic_biomass_change_rate']['#text']))
            nuclear_biomass_change_rate.append(float(volume_datum['nuclear_biomass_change_rate']['#text']))
            calcification_rate.append(float(volume_datum['calcification_rate']['#text']))

            fluid_fraction.append(float(volume_datum['fluid_fraction']['#text']))
            nuclear.append(float(volume_datum['nuclear']['#text']))
            calcified_fraction.append(float(volume_datum['calcified_fraction']['#text']))
            rel_rupture.append(None)
            total.append(float(volume_datum['total']['#text']))

    return phase_durations, fluid_change_rate, cytoplasmic_biomass_change_rate, nuclear_biomass_change_rate, \
           calcification_rate, fluid_fraction, nuclear, calcified_fraction, rel_rupture, total


def get_cycle_phenotypes(phenotypes, subdict, ppc):
    """
    The get_cycle_phenotypes function takes in three arguments: phenotypes, a dictionary of phenotype information
    already extracted,
    subdict, a sub-dictionary of the PhysiCell XML file that contains information about the current cell phenotype,
    and ppc, a dictionary of PhenoCellPy phenotype codes. This function extracts and processes cycle-specific phenotype
    data from the subdict and updates the phenotypes dictionary with the extracted data. If the cycle phenotype code is
    not in the ppc dictionary, the function sets the phenotype to the default "Simple Live" phenotype and issues a
    warning.

    The function extracts the duration and rate information for each cycle phase and stores it in the phenotypes
    dictionary, along with information about volume, fluid fraction, and calcification rate, if available. If volume
    information is not available, the function calculates phase durations based on rate information, or sets default
    values if rate information is not available.

    The function returns the updated phenotypes dictionary. If an error occurs during the extraction process, a
    ValueError is raised.

    :param phenotypes: dictionary with information about the phenotype models
    :param subdict: pcdict['cell_definitions']['cell_definition']
    :param ppc: codes of phenotypes
    :return: updated phenotypes dictionary
    """
    if subdict['phenotype']['cycle']['@code'] not in ppc.keys():
        message = f"WARNING: PhysiCell phenotype of code {subdict['phenotype']['cycle']['@code']}\n" \
                  f"not among PhenoCellPy's phenotypes. Falling back on Simple Live phenotype"
        warnings.warn(message)
        phenotypes[ppc["5"]] = None
        return phenotypes
    else:
        code_name = subdict['phenotype']['cycle']['@code']
        phenotype = ppc[code_name]
        if 'phase_transition_rates' in subdict['phenotype']['cycle'].keys():
            using_rates = True
            pheno_data = subdict['phenotype']['cycle']['phase_transition_rates']
        elif 'phase_durations' in subdict['phenotype']['cycle'].keys():
            using_rates = False
            pheno_data = subdict['phenotype']['cycle']['phase_durations']
        else:
            raise ValueError(f"Couldn't find phenotype phase transition data for "
                             f"{subdict['phenotype']['cycle']['name']}.\nIs this PhisiCell model valid?")

        rate_data = pheno_data['rate']
        if 'volume' in subdict['phenotype'].keys():
            volume_datum = subdict['phenotype']['volume']
            phase_durations, fluid_change_rate, cytoplasmic_biomass_change_rate, nuclear_biomass_change_rate, \
            calcification_rate, fluid_fraction, nuclear, calcified_fraction, rel_rupture, total = \
                get_cycle_rate_data(rate_data, volume_datum, using_rates)

            phenotypes[phenotype] = {"rate units": pheno_data['@units'],
                                     "phase durations": phase_durations,
                                     "fluid fraction": fluid_fraction,
                                     "fluid change rate": fluid_change_rate,
                                     "nuclear volume": nuclear,
                                     "cytoplasm biomass change rate": cytoplasmic_biomass_change_rate,
                                     "nuclear biomass change rate": nuclear_biomass_change_rate,
                                     "calcified fraction": calcified_fraction,
                                     "calcification rate": calcification_rate,
                                     "relative rupture volume": rel_rupture,
                                     "total": total}

        else:
            if 'phase_transition_rates' in subdict['phenotype']['cycle'].keys():
                using_rates = True
                pheno_data = subdict['phenotype']['cycle']['phase_transition_rates']
            elif 'phase_durations' in subdict['phenotype']['cycle'].keys():
                using_rates = False
                pheno_data = subdict['phenotype']['cycle']['phase_durations']
            if using_rates:
                phase_duration = 1 / float(rate_data['#text']) if float(rate_data['#text']) \
                    else 9e99
            else:
                phase_duration = float(rate_data['#text']) if float(rate_data['#text']) \
                    else 9e99
            fixed_duration = rate_data['@fixed_duration'].upper()
            duration_data = (fixed_duration, phase_duration)
            phase_durations = [duration_data]
            fluid_change_rate = [None]
            cytoplasmic_biomass_change_rate = [None]
            nuclear_biomass_change_rate = [None]
            calcification_rate = [None]
            fluid_fraction = [None]
            nuclear = [None]
            calcified_fraction = [None]
            phenotypes[phenotype] = {"rate units": pheno_data['@units'],
                                     "phase durations": phase_durations,
                                     "fluid fraction": fluid_fraction,
                                     "fluid change rate": fluid_change_rate,
                                     "nuclear volume": nuclear,
                                     "cytoplasm biomass change rate": cytoplasmic_biomass_change_rate,
                                     "nuclear biomass change rate": nuclear_biomass_change_rate,
                                     "calcified fraction": calcified_fraction,
                                     "calcification rate": calcification_rate,
                                     "relative rupture volume": [None]}
    return phenotypes


def get_death_phenotypes(phenotypes, subdict, ppc):
    """
    Extracts information about cell death phenotypes from a PhysiCell configuration subdictionary.

    The get_death_phenotypes function extracts information about cell death phenotypes from a PhysiCell configuration
    file, and returns a dictionary containing the updated cell phenotypes, including information about cell death. It
    takes in three arguments:

    * phenotypes: A dictionary containing information about cell phenotypes.
    * subdict: A sub-dictionary of the PhysiCell configuration file that contains information about the death models.
    * ppc: A dictionary that maps PhysiCell phenotype codes to PhenoCellPy phenotype names.

    The function returns a dictionary containing the updated cell phenotypes, including information about cell death.
    The function does not raise any exceptions.

    Parameters:
    -----------
        phenotypes : dict
            A dictionary containing information about cell phenotypes.
        subdict : dict
            A sub-dictionary of the PhysiCell configuration file that contains information about the death models.
        ppc : dict
            A dictionary that maps PhysiCell phenotype codes to PhenoCellPy phenotype names.

    Returns:
    -----------
        phenotypes : dict
            A dictionary containing the updated cell phenotypes, including information about cell death.

    """
    death_models = subdict['phenotype']['death']['model']
    if type(death_models) == list:
        for model in death_models:
            code_name = model['@code']
            if code_name not in ppc.keys():
                message = f"WARNING: PhysiCell phenotype of code {subdict['phenotype']['cycle']['@code']}\n" \
                          f"not among PhenoCellPy's phenotypes. Falling back on Standard apoptosis model phenotype"
                warnings.warn(message)
                phenotypes[ppc["100"]] = None
            else:
                phenotype = ppc[code_name]
                phenotypes[phenotype] = {"rate units": model['death_rate']['@units']}
                phase_durations = model['phase_durations']['duration']
                duration_data = []
                if type(phase_durations) == list:
                    for phasedur in phase_durations:
                        fixed = phasedur['@fixed_duration'].upper()
                        duration = float(phasedur['#text']) if float(phasedur['#text']) else 9e99
                        duration_data.append((fixed, duration))
                else:
                    fixed = phase_durations['@fixed_duration'].upper()
                    duration = float(phase_durations['#text']) if float(phase_durations['#text']) else 9e99
                    duration_data.append((fixed, duration))
                phenotypes[phenotype]["phase durations"] = duration_data

                biomass_chage_rates = model['parameters']
                if code_name == "100":  # apoptosis
                    phenotypes[phenotype]["fluid change rate"] = \
                        [float(biomass_chage_rates['unlysed_fluid_change_rate']['#text'])]
                    phenotypes[phenotype]["cytoplasm biomass change rate"] = \
                        [float(biomass_chage_rates['cytoplasmic_biomass_change_rate']['#text'])]
                    phenotypes[phenotype]["nuclear biomass change rate"] = \
                        [float(biomass_chage_rates['nuclear_biomass_change_rate']['#text'])]
                    phenotypes[phenotype]["calcification rate"] = \
                        [float(biomass_chage_rates['calcification_rate']['#text'])]
                    phenotypes[phenotype]["relative rupture volume"] = [None]
                elif code_name == "101":  # necrosis
                    phenotypes[phenotype]["fluid change rate"] = \
                        [float(biomass_chage_rates['unlysed_fluid_change_rate']['#text']),
                         float(biomass_chage_rates['lysed_fluid_change_rate']['#text'])]
                    phenotypes[phenotype]["cytoplasm biomass change rate"] = \
                        [float(biomass_chage_rates['cytoplasmic_biomass_change_rate']['#text']),
                         float(biomass_chage_rates['cytoplasmic_biomass_change_rate']['#text'])]
                    phenotypes[phenotype]["nuclear biomass change rate"] = \
                        [float(biomass_chage_rates['nuclear_biomass_change_rate']['#text']),
                         float(biomass_chage_rates['nuclear_biomass_change_rate']['#text'])]
                    phenotypes[phenotype]["calcification rate"] = \
                        [float(biomass_chage_rates['calcification_rate']['#text']),
                         float(biomass_chage_rates['calcification_rate']['#text'])]
                    phenotypes[phenotype]["relative rupture volume"] = [None, 2]
    else:
        model = death_models
        code_name = model['@code']
        if code_name not in ppc.keys():
            message = f"WARNING: PhysiCell phenotype of code {subdict['phenotype']['cycle']['@code']}\n" \
                      f"not among PhenoCellPy's phenotypes. Falling back on Standard apoptosis model phenotype"
            warnings.warn(message)
            phenotypes[ppc["100"]] = None
            return phenotypes
        else:
            phenotype = ppc[code_name]
            phenotypes[phenotype] = {"rate units": model['death_rate']['@units']}
            if 'phase_durations' in model.keys():
                phase_durations = model['phase_durations']['duration']
                duration_data = []
                if type(phase_durations) == list:
                    for phasedur in phase_durations:
                        fixed = phasedur['@fixed_duration'].upper()
                        duration = float(phasedur['#text'])
                        duration_data.append((fixed, duration))
                else:
                    fixed = phase_durations['@fixed_duration'].upper()
                    duration = float(phase_durations['#text'])
                    duration_data.append((fixed, duration))
                phenotypes[phenotype]["phase durations"] = duration_data
            else:
                phenotypes[phenotype]["phase durations"] = [(None, None)] * len(phenotypes[phenotype]["rate units"])

            if 'parameters' in model.keys():
                biomass_chage_rates = model['parameters']
                if code_name == "100":  # apoptosis
                    phenotypes[phenotype]["fluid change rate"] = \
                        [float(biomass_chage_rates['unlysed_fluid_change_rate']['#text'])]
                    phenotypes[phenotype]["cytoplasm biomass change rate"] = \
                        [float(biomass_chage_rates['cytoplasmic_biomass_change_rate']['#text'])]
                    phenotypes[phenotype]["nuclear biomass change rate"] = \
                        [float(biomass_chage_rates['nuclear_biomass_change_rate']['#text'])]
                    phenotypes[phenotype]["calcification rate"] = \
                        [float(biomass_chage_rates['calcification_rate']['#text'])]
                    phenotypes[phenotype]["relative rupture volume"] = [None]
                elif code_name == "101":  # necrosis
                    phenotypes[phenotype]["fluid change rate"] = \
                        [float(biomass_chage_rates['unlysed_fluid_change_rate']['#text']),
                         float(biomass_chage_rates['lysed_fluid_change_rate']['#text'])]
                    phenotypes[phenotype]["cytoplasm biomass change rate"] = \
                        [float(biomass_chage_rates['cytoplasmic_biomass_change_rate']['#text']),
                         float(biomass_chage_rates['cytoplasmic_biomass_change_rate']['#text'])]
                    phenotypes[phenotype]["nuclear biomass change rate"] = \
                        [float(biomass_chage_rates['nuclear_biomass_change_rate']['#text']),
                         float(biomass_chage_rates['nuclear_biomass_change_rate']['#text'])]
                    phenotypes[phenotype]["calcification rate"] = \
                        [float(biomass_chage_rates['calcification_rate']['#text']),
                         float(biomass_chage_rates['calcification_rate']['#text'])]
                    phenotypes[phenotype]["relative rupture volume"] = [None, 2]
    return phenotypes


def get_cell_phenotypes(subdict, ppc=_physicell_phenotype_codes):
    """
    Extracts the cell phenotypes for a given cell from a pcdict['cell_definitions']['cell_definition'] subdictionary.

    Parameters:
    -----------
    subdict : dict
        A dictionary containing information about the cell.
    ppc : dict, optional
        A dictionary of PhysiCell phenotype codes, where the keys are the codes and the values are the corresponding
        phenotype names. The default is _physicell_phenotype_codes.

    Returns:
    --------
    phenotypes : dict or None
        A dictionary containing the cell phenotypes and their values. Returns None if the given subdictionary does not
        contain any phenotypes.
    pheno_names : list or None
        A list of the names of the cell phenotypes. Returns None if the given subdictionary does not contain any
        phenotypes.

    Notes:
    ------
    The phenotypes are extracted based on the keys 'cycle' and 'death' in the subdictionary. If a key is present, the
    corresponding phenotype values are extracted using the helper functions `get_cycle_phenotypes()` and
    `get_death_phenotypes()`.

    """
    phenotypes = {}
    if 'phenotype' not in subdict.keys():
        return None, None
    if 'cycle' in subdict['phenotype'].keys():
        phenotypes = get_cycle_phenotypes(phenotypes, subdict, ppc)
    if 'death' in subdict['phenotype'].keys():
        phenotypes = get_death_phenotypes(phenotypes, subdict, ppc)
    pheno_names = list(phenotypes.keys())
    return phenotypes, pheno_names


def get_custom_data(subdict):
    """
    Extracts the custom data for a given cell from a pcdict['cell_definitions']['cell_definition'] subdictionary.

    Parameters:
    -----------
    subdict : dict
        A dictionary containing information about the cell.

    Returns:
    --------
    custom_data : dict or None
        A dictionary containing the custom data for the given cell. Returns None if the given subdictionary does not
        contain custom data.
    """
    if "custom_data" in subdict.keys():
        return subdict["custom_data"]

    return None


def get_cell_constraints(pcdict, space_unit, minimum_volume=8):
    """
    Extracts cell constraints from the given PhysiCell pcdict.

    Parameters:
    -----------
    pcdict : dict
        Dictionary created from parsing PhysiCell XML. Must contain a "cell_definitions" key that maps
        to a dictionary with a "cell_definition" key. This key should contain a list of dictionaries, each of which
        represents a Cell Type.
    space_unit : float
        A scaling factor for the simulation's spatial units. All volumes extracted from pcdict will be multiplied by
        this factor raised to the power of the dimensionality of the simulation space.
    minimum_volume : float, optional
        The minimum volume allowed for any cell in pixels. If a cell's volume falls below this threshold after scaling,
        the translator will reconvert space so that the minimum cell volume is  equal to this threshold. Defaults to 8.

    Returns:
    --------
    constraints : dict
        A dictionary containing the constraints for each Cell Type found in pcdict. Each key is a Cell Type name
        (converted to an underscore-delimited string), and each value is a dictionary containing information about
        that Cell Type's volume, mechanics, custom data, and phenotypes.
    any_below : bool
        A boolean indicating whether any cells had volumes that fell below minimum_volume after scaling.
    volumes : list
        A list containing the scaled volumes of each Cell Type found in pcdict.
    minimum_volume : float
        The minimum volume allowed for any cell, after scaling.

    Raises:
    -------
    UserWarning
        If a Cell Type's volume is missing a unit or value in pcdict, or if the scaled volume falls below
        minimum_volume.
    """

    constraints = {}
    any_below = False
    volumes = []
    for child in pcdict['cell_definitions']['cell_definition']:
        ctype = child['@name'].replace(" ", "_")
        constraints[ctype] = {}
        volume, unit = get_cell_volume(child)
        if volume is None or unit is None:
            message = f"WARNING: cell volume for cell type {ctype} either doesn't have a unit \n(unit found: {unit}) " \
                      f"or" \
                      f" doesn't have a value (value found: {volume}). \nSetting the volume to be the minimum volume, " \
                      f"{minimum_volume}"
            warnings.warn(message)
            volume = None
            unit = "not specified"
            dim = 3
        else:
            dim = int(unit.split("^")[-1])
        if volume is None:
            volumepx = None
        else:
            volumepx = volume * (space_unit ** dim)
        below, minimum_volume = check_below_minimum_volume(volumepx, minimum=minimum_volume)
        constraints[ctype]["volume"] = {f"volume ({unit})": volume,
                                        "volume (pixels)": volumepx}
        volumes.append(volumepx)
        if below:
            any_below = True
            message = f"WARNING: converted cell volume for cell type {ctype} is below {minimum_volume}. Converted volume " \
                      f"{volumepx}. \nIf cells are too small in CC3D they do not behave in a biological manner and may " \
                      f"disapear. \nThis program will enforce that: 1) the volume proportions stay as before; 2) the " \
                      f"lowest cell volume is {minimum_volume}"
            warnings.warn(message)
        constraints[ctype]["mechanics"] = get_cell_mechanics(child)
        constraints[ctype]["custom_data"] = get_custom_data(child)
        constraints[ctype]["phenotypes"], constraints[ctype]["phenotypes_names"] = get_cell_phenotypes(child)

    return constraints, any_below, volumes, minimum_volume


def get_space_time_from_diffusion(unit):
    """
    This function takes a unit of measurement string in the format of 'spaceunit/timeunit', where spaceunit may contain
    an exponent, and returns a tuple with the extracted space unit and time unit as separate strings.

    Parameters:
    ----------
        unit (string): The unit of measurement in the format 'spaceunit/timeunit', where spaceunit may contain an
        exponent.

    Returns:
    -------
        A tuple containing two strings: the extracted space unit (without the exponent) and the time unit.
    """
    parts = unit.split("/")
    timeunit = parts[-1]
    spaceunit = parts[0].split("^")[0]
    return spaceunit, timeunit


def get_secretion_uptake(pcdict):
    """
    Extracts secretion data from the input pcdict (PhysiCell XML parsed into a dictionary) and returns the extracted
    data as a dictionary.

    get_secretion takes a single argument, pcdict, which is a Python dictionary created from converting a PhysiCell XML
    into a dictionary.

    The function initializes an empty dictionary sec_up_data, which will store the parsed secretion data. The code loops
    through the children of pcdict['cell_definitions']['cell_definition'] and checks if the child has the key
    'secretion' in its 'phenotype' dictionary. If not, the loop goes to the next child.

    For each child that has the 'secretion' key, the code extracts the child cell type as ctype. It then initializes a
    dictionary sec_up_data[ctype] to store secretion data for this cell type.

    The code then extracts a list of substrates and their secretion data for the given cell type and diffusing element.
    It stores this data in `sec_list` (`sec_list = child['phenotype']['secretion']['substrate']`).
    The code handles two cases: either `sec_list` is a list of multiple substrates (diffusing elements), or it is a
    single dictionary representing a single substrate.

    For each substrate, the code extracts its name, and then extracts various secretion data for that
    substrate, such as secretion_rate, uptake_rate, and net_export, if they exist. These values are all converted to
    floats if they exist, and default to 0 if they do not. The code also extracts the units for each of these values,
    which are stored as strings.

    All of this data is then stored in the sec_up_data dictionary for the given child type and substrate.

    Once all children with 'secretion' keys have been processed, the sec_up_data dictionary is returned.

    Parameters
    ----------
    pcdict : dict
        A dictionary representing a PhysiCell simulation.

    Returns
    -------
    dict
        A dictionary containing secretion data for each cell type and substrate.
    """

    sec_up_data = {}
    for child in pcdict['cell_definitions']['cell_definition']:
        if 'secretion' not in child['phenotype'].keys():
            continue
        ctype = child['@name'].replace(" ", "_")
        sec_up_data[ctype] = {}
        sec_list = child['phenotype']['secretion']['substrate']
        if type(sec_list) == list:
            for sec in sec_list:
                substrate = sec["@name"].replace(" ", "_")
                sec_up_data[ctype][substrate] = {}
                sec_up_data[ctype][substrate]['secretion_rate'] = float(
                    sec['secretion_rate']['#text']) if 'secretion_rate' in sec.keys() else 0
                sec_up_data[ctype][substrate]['secretion_unit'] = sec['secretion_rate'][
                    '@units'] if 'secretion_rate' in sec.keys() else "None"
                sec_up_data[ctype][substrate]['secretion_target'] = float(
                    sec['secretion_target']['#text']) if 'secretion_target' in sec.keys() else 0
                sec_up_data[ctype][substrate]['uptake_rate'] = float(
                    sec['uptake_rate']['#text']) if 'uptake_rate' in sec.keys() else 0
                sec_up_data[ctype][substrate]['uptake_unit'] = sec['uptake_rate'][
                    '@units'] if 'uptake_rate' in sec.keys() else "None"
                sec_up_data[ctype][substrate]['net_export'] = float(
                    sec['net_export_rate']['#text']) if 'net_export_rate' in sec.keys() else 0
                sec_up_data[ctype][substrate]['net_export_unit'] = sec['net_export_rate'][
                    '@units'] if 'net_export_rate' in sec.keys() else "None"
        else:
            sec = sec_list
            substrate = sec["@name"].replace(" ", "_")
            sec_up_data[ctype][substrate] = {}
            sec_up_data[ctype][substrate]['secretion_rate'] = float(
                sec['secretion_rate']['#text']) if 'secretion_rate' in sec.keys() else 0
            sec_up_data[ctype][substrate]['secretion_unit'] = sec['secretion_rate'][
                '@units'] if 'secretion_rate' in sec.keys() else "None"
            sec_up_data[ctype][substrate]['secretion_target'] = float(
                sec['secretion_target']['#text']) if 'secretion_target' in sec.keys() else 0
            sec_up_data[ctype][substrate]['uptake_rate'] = float(
                sec['uptake_rate']['#text']) if 'uptake_rate' in sec.keys() else 0
            sec_up_data[ctype][substrate]['uptake_unit'] = sec['uptake_rate'][
                '@units'] if 'uptake_rate' in sec.keys() else "None"
            sec_up_data[ctype][substrate]['net_export'] = float(
                sec['net_export_rate']['#text']) if 'net_export_rate' in sec.keys() else 0
            sec_up_data[ctype][substrate]['net_export_unit'] = sec['net_export_rate'][
                '@units'] if 'net_export_rate' in sec.keys() else "None"
    return sec_up_data


def get_microenvironment(pcdict, space_factor, space_unit, time_factor, time_unit, autoconvert_time=True,
                         autoconvert_space=True, space_convs=_space_convs, time_convs=_time_convs,
                         steady_state_threshold=1000):
    """
    Gets data from diffusing elements defined in PhysiCell and converts it into CompuCell3D ready values

    Given a dictionary with the simulation setup `pcdict`, the `space_factor` and `time_factor` to use for unit
    conversion, the desired `space_unit` and `time_unit` in which to express the diffusion and decay coefficients
    respectively, the function retrieves the microenvironment from `pcdict` and returns a dictionary with the diffusion
    and decay coefficient values converted to CC3D units.

    If automatic unit conversion is not possible, a warning message will be printed to the console and automatic unit
    conversion will be disabled for that variable.

    Parameters:
    -----------
        pcdict : dict
            The dictionary containing the PhysiCell simulation setup.
        space_factor : float
            The factor to use for space unit conversion.
        space_unit : str
            The unit used in PhysiCell for space
        time_factor : float
            The factor to use for time unit conversion.
        time_unit : str
            The unit used in PhysiCell for time
        autoconvert_time : bool, optional
            Whether to automatically convert time units into CC3D time units. Default is True.
        autoconvert_space : bool, optional
            Whether to automatically convert space units into CC3D space units. Default is True
        space_convs : dict, optional
            A dictionary of known space unit conversions. Default is _space_convs.
        time_convs : dict, optional
            A dictionary of known time unit conversions. Default is _time_convs.
        steady_state_threshold : float, optional
            The steady state threshold value above which the translator will set the diffusion solver  for that chemical
            to be the steady state solver. Default is 1000.

    Returns:
        dict: A dictionary containing the diffusion and decay coefficients, initial condition, and boundary
            conditions for each diffusing element in the microenvironment.

    Raises:
        ValueError: If the given space unit or time unit is not recognized or not convertible.
    """
    diffusing_elements = {}
    fields = pcdict['microenvironment_setup']['variable']
    for subel in fields:

        auto_s_this = autoconvert_space
        auto_t_this = autoconvert_time

        diffusing_elements[subel['@name']] = {}
        diffusing_elements[subel['@name']]["use_steady_state"] = False
        diffusing_elements[subel['@name']]["concentration_units"] = subel["@units"]
        diffusing_elements[subel['@name']]["D_w_units"] = \
            float(subel['physical_parameter_set']['diffusion_coefficient']['#text'])

        this_space, this_time = get_space_time_from_diffusion(subel['physical_parameter_set']['diffusion_coefficient']
                                                              ['@units'])

        if this_space != space_unit:
            message = f"WARNING: space unit found in diffusion coefficient of {subel['@name']} does not match" \
                      f"space unit found while converting <overall>:\n\t<overall>:{space_unit};\n\t{subel['@name']}:" \
                      f"{this_space}" \
                      f"\nautomatic space-unit conversion for {subel['@name']} disabled"
            warnings.warn(message)
            auto_s_this = False
            # space_conv_factor = 1
        if this_time != time_unit:
            message = f"WARNING: time unit found in diffusion coefficient of {subel['@name']} does not match" \
                      f"space unit found while converting <overall>:\n\t<overall>:{time_unit};" \
                      f"\n\t{subel['@name']}:{this_time}" \
                      f"\nautomatic time-unit conversion for {subel['@name']} disabled"
            warnings.warn(message)
            auto_t_this = False
            # time_conv_factor = 1

        diffusing_elements[subel['@name']]["auto"] = (auto_s_this, auto_t_this)

        if auto_s_this:
            space_conv_factor = space_factor
        else:
            space_conv_factor = 1

        if auto_t_this:
            time_conv_factor = time_factor
        else:
            time_conv_factor = 1

        D = diffusing_elements[subel['@name']]["D_w_units"] * space_conv_factor * space_conv_factor / time_conv_factor
        diffusing_elements[subel['@name']]["D"] = D
        if D > steady_state_threshold:
            diffusing_elements[subel['@name']]["use_steady_state"] = True

        # [cc3dds] = pixel/unit
        # [cc3dds] * unit = pixel -> pixel^2 = ([cc3dds] * unit)^2
        # [cc3ddt] = MCS/unit
        # [cc3ddt] * unit = MCS -> 1/MCS = 1/([cc3ddt] * unit)

        if auto_s_this and auto_t_this:
            diffusing_elements[subel['@name']]["D_conv_factor_text"] = f"1 pixel^2/MCS" \
                                                                       f" = {space_conv_factor * space_conv_factor / time_conv_factor}" \
                                                                       f" {space_unit}^2/{time_unit}"
            diffusing_elements[subel['@name']]["D_conv_factor"] = space_conv_factor * space_conv_factor / \
                                                                  time_conv_factor
            diffusing_elements[subel['@name']]["D_og_unit"] = f"{space_unit}^2/{time_unit}"
        else:
            diffusing_elements[subel['@name']]["D_conv_factor_text"] = "disabled autoconversion"
            diffusing_elements[subel['@name']]["D_conv_factor"] = 1
            diffusing_elements[subel['@name']]["D_og_unit"] = "disabled autoconversion, not known"
        diffusing_elements[subel['@name']]["gamma_w_units"] = float(subel['physical_parameter_set']['decay_rate']
                                                                    ['#text'])
        gamma = diffusing_elements[subel['@name']]["gamma_w_units"] / time_conv_factor
        diffusing_elements[subel['@name']]["gamma"] = gamma
        if auto_t_this:
            diffusing_elements[subel['@name']][
                "gamma_conv_factor_text"] = f"1/MCS = {1 / time_conv_factor} 1/{time_unit}"
            diffusing_elements[subel['@name']]["gamma_conv_factor"] = 1 / time_conv_factor
            diffusing_elements[subel['@name']]["gamma_og_unit"] = f" 1/{time_unit}"
        else:
            diffusing_elements[subel['@name']]["gamma_conv_factor_text"] = "disabled autoconversion"
            diffusing_elements[subel['@name']]["gamma_conv_factor"] = 1
            diffusing_elements[subel['@name']]["gamma_og_unit"] = "disabled autoconversion, not known"

        diffusing_elements[subel['@name']]["initial_condition"] = subel['initial_condition']['#text']

        diffusing_elements[subel['@name']]["dirichlet"] = subel['Dirichlet_boundary_condition']['@enabled']
        diffusing_elements[subel['@name']]["dirichlet_value"] = float(subel['Dirichlet_boundary_condition']['#text'])

    return diffusing_elements
