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
    secretion_comment = ''
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


def convert_secretion_data(sec_dict, time_conv, pctimeunit):
    if not sec_dict:
        return {}

    new_sec_dict = sec_dict

    secretion_comment = '#WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not. \n#The ' \
                        'translating program attempts to implement it, but it may not be a 1 to 1 conversion.'
    uptake_comment = '#WARNING: To avoid negative concentrations, in CompuCell3D uptake is "bounded." \n# If the amount' \
                     ' that would be uptaken is larger than the value at that pixel,\n# the uptake will be a set ratio ' \
                     'of the amount available.\n# The conversion program uses 1 as the ratio,\n# you may want to ' \
                     'revisit this.'
    for ctype in sec_dict.keys():
        type_sec = sec_dict[ctype]
        new_type_sec = type_sec
        for field, data in type_sec.items():
            unit = data['secretion_unit']

            mcs_secretion_rate, extra_sec_comment = convert_secretion_rate(data['secretion_rate'], unit, time_conv,
                                                                           pctimeunit)
            data['secretion_rate_MCS'] = mcs_secretion_rate
            data['secretion_comment'] = secretion_comment + extra_sec_comment

            unit = data['net_export_unit']

            mcs_net_secretion_rate, extra_net_sec_comment = convert_net_secretion(data['net_export'], unit, time_conv,
                                                                                  pctimeunit)
            data['net_export_MCS'] = mcs_net_secretion_rate
            data['net_secretion_comment'] = extra_net_sec_comment

            mcs_uptake_rate, extra_up_comment = convert_uptake_rate(data['uptake_rate'], data['uptake_unit'], time_conv,
                                                                    pctimeunit)

            data['uptake_rate_MCS'] = mcs_uptake_rate
            data['uptake_comment'] = uptake_comment + extra_up_comment

            new_type_sec[field] = data
        new_sec_dict[ctype] = new_type_sec
    return new_sec_dict