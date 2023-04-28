import warnings

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


def convert_secretion_rate(rate, unit, time_conv, pctimeunit, time_convs=_time_convs):
    """
    Converts a secretion rate from its original units to MCS units.

    If the original units cannot be determined, no conversion is performed and a warning message is displayed.
    If the original units are the same as the main PhysiCell time unit, the conversion factor is simply the
    simulation time conversion factor.
    If the original units are different from the main PhysiCell time unit, a conversion factor is determined from
    the 'time_convs' dictionary and used to convert the secretion rate to minutes. The converted rate is then
    converted to the main PhysiCell time unit, and finally to MCS units using the simulation time conversion factor.

    Parameters
    ----------
    rate : float
       The secretion rate value to be converted.
    unit : str
       The original units of the secretion rate value.
    time_conv : float
       Conversion factor from original simulation time units to minutes.
    pctimeunit : str
       The main PhysiCell time unit used in the simulation.
    time_convs : dict, optional
       Dictionary of conversion factors from other time units to minutes, by default _time_convs.

    Returns
    -------
    float
       The converted secretion rate value in MCS units.
    str
       A warning message with additional information about the conversion process.
    """
    secretion_comment = ''
    if unit is None:
        mcs_rate = rate
        message = f"WARNING: couldn't find units for this rate!\nAutomatic conversion" \
                  f" of " \
                  f"this rate is disabled."
        secretion_comment += "\n#" + message.replace("\n", "\n#")
        warnings.warn(message)
        return mcs_rate, secretion_comment
    if pctimeunit in unit:  # if it's the same as the "main" time unit
        mcs_rate = rate / time_conv
        return mcs_rate, secretion_comment
    else:
        tu = unit.split('/')[-1]
        if tu not in time_convs.keys():
            message = f"WARNING: Secretion 1/(rate unit) = {tu} not found in {time_convs.keys()}.\nAutomatic conversion" \
                      f" of " \
                      f"this rate is disabled."
            secretion_comment += "\n#" + message.replace("\n", "\n#")
            warnings.warn(message)
            mcs_rate = rate
            return mcs_rate, secretion_comment
        else:
            message = f"WARNING: Secretion 1/(rate unit) = {tu} is not the main PhysiCell time unit {pctimeunit}." \
                      f"\nTherefore, the automatic conversion may be incorrect."
            secretion_comment += "\n#" + message.replace("\n", "\n#")
            warnings.warn(message)
            rate_minutes = rate / time_convs[tu]
            rate_pctime = rate_minutes / time_conv[pctimeunit]
            mcs_rate = rate_pctime / time_conv
            return mcs_rate, secretion_comment


def convert_uptake_rate(rate, unit, time_conv, pctimeunit, time_convs=_time_convs):
    uptake_comment = ''
    if pctimeunit in unit:
        mcs_rate = rate / time_conv
        return mcs_rate, uptake_comment
    else:
        tu = unit.split('/')[-1]
        if tu not in time_convs.keys():
            message = f"WARNING: Uptake 1/(rate unit) = {tu} not found in {time_convs.keys()}.\nAutomatic " \
                      f"conversion of this rate is disabled."
            warnings.warn(message)
            uptake_comment += "#" + message.replace("\n", "\n#")
            mcs_rate = rate
            return mcs_rate, uptake_comment
        else:
            message = f"WARNING: Uptake 1/(rate unit) = {tu} is not the main PhysiCell time unit {pctimeunit}." \
                      f"\nTherefore, the automatic conversion may be incorrect"
            warnings.warn(message)
            uptake_comment += "#" + message.replace("\n", "\n#")

            rate_minutes = rate / time_convs[tu]
            rate_pctime = rate_minutes / time_conv[pctimeunit]
            mcs_rate = rate_pctime / time_conv
            return mcs_rate, uptake_comment


def convert_net_secretion(rate, unit, time_conv, pctimeunit, time_convs=_time_convs):
    """
    Convert a net secretion rate from a given time unit to 1/MCS.

    If the unit of the rate is not expressed in `pctimeunit`, the function attempts to convert it automatically using
    the `time_convs` dictionary. If the time unit of the rate is not found in `time_convs` a warning message is
    generated and no conversion is performed. If the time unit is found but
    is not the same as the main PhysiCell time unit, a warning message is generated and conversion is performed.
    The warning message is returned as part of the `net_comment` string.

    Parameters
    ----------
    rate : float
        The net secretion rate in the given time unit.
    unit : str
        The unit of the rate, expressed as a string (e.g., "1/min", "1/hour", etc.).
    time_conv : float
        The conversion factor from the given time unit to 1/MCS
    pctimeunit : str
        The main PhysiCell time unit, expressed as a string (e.g., "min").
    time_convs : dict, optional
        A dictionary containing conversion factors for different time units. The keys are time units expressed as
        strings (e.g., "min", "hour", etc.), and the values are the corresponding conversion factors to minutes.

    Returns
    -------
    mcs_rate : float
        The net secretion rate converted to 1/MCS
    net_comment : str
        A warning message (if any) about the conversion, expressed as a string with comment marks ('#') at the
        beginning of each line.
    """
    net_comment = ''
    if pctimeunit in unit:
        mcs_rate = rate / time_conv
        return mcs_rate, net_comment
    else:
        tu = unit.split('/')[-1]
        if tu not in time_convs.keys():
            message = f"WARNING: Time component of net secretion unit ({unit}) not found in {time_convs.keys()}." \
                      f"\nAutomatic conversion of this rate is disabled."
            warnings.warn(message)
            net_comment += "#" + message.replace("\n", "\n#")
            mcs_rate = rate
            return mcs_rate, net_comment
        else:
            message = f"WARNING: Time component of net secretion unit ({unit}) is not the main PhysiCell time unit " \
                      f"{pctimeunit}. \nTherefore, the automatic conversion may be incorrect"
            warnings.warn(message)
            net_comment += "#" + message.replace("\n", "\n#")
            rate_minutes = rate / time_convs[tu]
            rate_pctime = rate_minutes / time_conv[pctimeunit]
            mcs_rate = rate_pctime / time_conv
            return mcs_rate, net_comment


def convert_secretion_uptake_data(sec_dict, time_conv, pctimeunit):
    """
    Convert secretion data from PhysiCell to CompuCell3D format.

    This function converts the secretion data parsed from PhysiCell to CompuCell3D python commands to perform
    secretion. It also adds comment to indicate any potential discrepancies between the two formats, such as
    differences in the handling of target secretion or uptake bounds. The resulting dictionary is returned.

    Parameters
    ----------
    sec_dict : dict
        Dictionary containing parsed secretion data from PhysiCell.
    time_conv : float
        Conversion factor for time units.
    pctimeunit : str
        PhysiCell time unit.

    Returns
    -------
    dict
        Dictionary containing CompuCell3D secretion python commands.

    Notes
    -----


    """
    if not sec_dict:
        return {}

    new_sec_dict = sec_dict

    secretion_comment = '#WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not. \n#The ' \
                        'translating program attempts to implement it, but it may not be a 1 to 1 conversion.'
    uptake_comment = '#WARNING: To avoid negative concentrations, in CompuCell3D uptake is "bounded." \n# If the amount' \
                     ' that would be uptaken is larger than the value at that pixel,\n# the uptake will be a set ratio ' \
                     'of the amount available.\n# The conversion program uses 1 as the ratio,\n# you may want to ' \
                     'revisit this.'

    # secretion in physicell is
    # secretion rate * (target amount - amount at cell) + net secretion
    # looking at units that is correct:
    # < secretion_rate units = "1/min" > 0 < / secretion_rate >
    # < secretion_target units = "substrate density" > 1 < / secretion_target >
    # < uptake_rate units = "1/min" > 0 < / uptake_rate >
    # < net_export_rate units = "total substrate/min" > 0 < / net_export_rate >

    for ctype in sec_dict.keys():
        type_sec = sec_dict[ctype]
        new_type_sec = type_sec
        for field, data in type_sec.items():
            unit = data['secretion_unit'] if 'secretion_unit' in data.keys() else None

            mcs_secretion_rate, extra_sec_comment = convert_secretion_rate(data['secretion_rate'], unit, time_conv,
                                                                           pctimeunit)
            data['secretion_rate_MCS'] = mcs_secretion_rate
            data['secretion_comment'] = secretion_comment + extra_sec_comment

            # data["secretion_target"] = get_secretion_target()

            unit = data['net_export_unit'] if 'net_export_unit' in data.keys() else None

            mcs_net_secretion_rate, extra_net_sec_comment = convert_net_secretion(data['net_export'], unit, time_conv,
                                                                                  pctimeunit)
            data['net_export_MCS'] = mcs_net_secretion_rate
            data['net_secretion_comment'] = extra_net_sec_comment

            mcs_uptake_rate, extra_up_comment = convert_uptake_rate(data['uptake_rate'], data['uptake_unit'], time_conv,
                                                                    pctimeunit)

            data['uptake_rate_MCS'] = mcs_uptake_rate
            data['uptake_comment'] = uptake_comment + extra_up_comment
            if "secretion_target" not in data.keys():
                data["secretion_target"] = 0
            new_type_sec[field] = data

        new_sec_dict[ctype] = new_type_sec

    return new_sec_dict
