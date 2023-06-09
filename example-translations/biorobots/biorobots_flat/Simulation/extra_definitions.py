cell_constraints = {
    'default': {
        'volume': {
            'volume (micron^3)': 2494.0,
            'volume (pixels)': 8.105500000000003},
        'mechanics': {
            'cell_cell_adhesion_strength': {
                'units': 'micron/min',
                'value': 0.4},
            'cell_cell_repulsion_strength': {
                'units': 'micron/min',
                'value': 10.0},
            'relative_maximum_adhesion_distance': {
                'units': 'dimensionless',
                'value': 1.25}},
        'custom_data': OrderedDict(
            [
                ('receptor',
                 OrderedDict(
                     [
                         ('@units',
                          'dimensionless'),
                         ('#text',
                          '0.0')]))]),
        'phenotypes': {
            'Simple Live': {
                'rate units': '1/min',
                'phase durations': [
                    ('FALSE',
                     9e+99)],
                'fluid fraction': [0.75],
                'fluid change rate': [0.05],
                'nuclear volume': [540.0],
                'cytoplasm biomass change rate': [0.0045],
                'nuclear biomass change rate': [0.0055],
                'calcified fraction': [0.0],
                'calcification rate': [0.0],
                'relative rupture volume': [None],
                'total': [2494.0]},
            'Standard apoptosis model': {
                'rate units': '1/min',
                'phase durations': [
                    ('TRUE',
                     516.0)],
                'fluid change rate': [0.05],
                'cytoplasm biomass change rate': [0.0166667],
                'nuclear biomass change rate': [0.00583333],
                'calcification rate': [0.0],
                'relative rupture volume': [None]},
            'Standard necrosis model': {
                'rate units': '1/min',
                'phase durations': [
                    ('TRUE',
                     9e+99),
                    ('TRUE',
                     86400.0)],
                'fluid change rate': [
                    0.05,
                    0.0],
                'cytoplasm biomass change rate': [
                    0.0166667,
                    0.0166667],
                'nuclear biomass change rate': [
                    0.00583333,
                    0.00583333],
                'calcification rate': [
                    0.0,
                    0.0],
                'relative rupture volume': [
                    None,
                    2]}},
        'phenotypes_names': [
            'Simple Live',
            'Standard apoptosis model',
            'Standard necrosis model']},
    'director_cell': {
        'volume': {
            'volume (micron^3)': 2494.0,
            'volume (pixels)': 8.105500000000003},
        'mechanics': {
            'cell_cell_adhesion_strength': {
                'units': 'micron/min',
                'value': 0.4},
            'cell_cell_repulsion_strength': {
                'units': 'micron/min',
                'value': 10.0},
            'relative_maximum_adhesion_distance': {
                'units': 'dimensionless',
                'value': 1.25}},
        'custom_data': OrderedDict(
            [
                ('receptor',
                 OrderedDict(
                     [
                         ('@units',
                          'dimensionless'),
                         ('#text',
                          '0.0')]))]),
        'phenotypes': {
            'Simple Live': {
                'rate units': '1/min',
                'phase durations': [
                    ('FALSE',
                     9e+99)],
                'fluid fraction': [0.75],
                'fluid change rate': [0.05],
                'nuclear volume': [540.0],
                'cytoplasm biomass change rate': [0.0045],
                'nuclear biomass change rate': [0.0055],
                'calcified fraction': [0.0],
                'calcification rate': [0.0],
                'relative rupture volume': [None],
                'total': [2494.0]},
            'Standard apoptosis model': {
                'rate units': '1/min',
                'phase durations': [
                    ('TRUE',
                     516.0)],
                'fluid change rate': [0.05],
                'cytoplasm biomass change rate': [0.0166667],
                'nuclear biomass change rate': [0.00583333],
                'calcification rate': [0.0],
                'relative rupture volume': [None]},
            'Standard necrosis model': {
                'rate units': '1/min',
                'phase durations': [
                    ('TRUE',
                     9e+99),
                    ('TRUE',
                     86400.0)],
                'fluid change rate': [
                    0.05,
                    0.0],
                'cytoplasm biomass change rate': [
                    0.0166667,
                    0.0166667],
                'nuclear biomass change rate': [
                    0.00583333,
                    0.00583333],
                'calcification rate': [
                    0.0,
                    0.0],
                'relative rupture volume': [
                    None,
                    2]}},
        'phenotypes_names': [
            'Simple Live',
            'Standard apoptosis model',
            'Standard necrosis model']},
    'cargo_cell': {
        'volume': {
            'volume (micron^3)': 2494.0,
            'volume (pixels)': 8.105500000000003},
        'mechanics': {
            'cell_cell_adhesion_strength': {
                'units': 'micron/min',
                'value': 0.4},
            'cell_cell_repulsion_strength': {
                'units': 'micron/min',
                'value': 10.0},
            'relative_maximum_adhesion_distance': {
                'units': 'dimensionless',
                'value': 1.25}},
        'custom_data': OrderedDict(
            [
                ('receptor',
                 OrderedDict(
                     [
                         ('@units',
                          'dimensionless'),
                         ('#text',
                          '1.0')]))]),
        'phenotypes': {
            'Simple Live': {
                'rate units': '1/min',
                'phase durations': [
                    ('FALSE',
                     9e+99)],
                'fluid fraction': [0.75],
                'fluid change rate': [0.05],
                'nuclear volume': [540.0],
                'cytoplasm biomass change rate': [0.0045],
                'nuclear biomass change rate': [0.0055],
                'calcified fraction': [0.0],
                'calcification rate': [0.0],
                'relative rupture volume': [None],
                'total': [2494.0]},
            'Standard apoptosis model': {
                'rate units': '1/min',
                'phase durations': [
                    ('TRUE',
                     516.0)],
                'fluid change rate': [0.05],
                'cytoplasm biomass change rate': [0.0166667],
                'nuclear biomass change rate': [0.00583333],
                'calcification rate': [0.0],
                'relative rupture volume': [None]},
            'Standard necrosis model': {
                'rate units': '1/min',
                'phase durations': [
                    ('TRUE',
                     9e+99),
                    ('TRUE',
                     86400.0)],
                'fluid change rate': [
                    0.05,
                    0.0],
                'cytoplasm biomass change rate': [
                    0.0166667,
                    0.0166667],
                'nuclear biomass change rate': [
                    0.00583333,
                    0.00583333],
                'calcification rate': [
                    0.0,
                    0.0],
                'relative rupture volume': [
                    None,
                    2]}},
        'phenotypes_names': [
            'Simple Live',
            'Standard apoptosis model',
            'Standard necrosis model']},
    'worker_cell': {
        'volume': {
            'volume (micron^3)': 2494.0,
            'volume (pixels)': 8.105500000000003},
        'mechanics': {
            'cell_cell_adhesion_strength': {
                'units': 'micron/min',
                'value': 0.4},
            'cell_cell_repulsion_strength': {
                'units': 'micron/min',
                'value': 10.0},
            'relative_maximum_adhesion_distance': {
                'units': 'dimensionless',
                'value': 1.25}},
        'custom_data': OrderedDict(
            [
                ('receptor',
                 OrderedDict(
                     [
                         ('@units',
                          'dimensionless'),
                         ('#text',
                          '1.0')]))]),
        'phenotypes': {
            'Simple Live': {
                'rate units': '1/min',
                'phase durations': [
                    ('FALSE',
                     9e+99)],
                'fluid fraction': [0.75],
                'fluid change rate': [0.05],
                'nuclear volume': [540.0],
                'cytoplasm biomass change rate': [0.0045],
                'nuclear biomass change rate': [0.0055],
                'calcified fraction': [0.0],
                'calcification rate': [0.0],
                'relative rupture volume': [None],
                'total': [2494.0]},
            'Standard apoptosis model': {
                'rate units': '1/min',
                'phase durations': [
                    ('TRUE',
                     516.0)],
                'fluid change rate': [0.05],
                'cytoplasm biomass change rate': [0.0166667],
                'nuclear biomass change rate': [0.00583333],
                'calcification rate': [0.0],
                'relative rupture volume': [None]},
            'Standard necrosis model': {
                'rate units': '1/min',
                'phase durations': [
                    ('TRUE',
                     9e+99),
                    ('TRUE',
                     86400.0)],
                'fluid change rate': [
                    0.05,
                    0.0],
                'cytoplasm biomass change rate': [
                    0.0166667,
                    0.0166667],
                'nuclear biomass change rate': [
                    0.00583333,
                    0.00583333],
                'calcification rate': [
                    0.0,
                    0.0],
                'relative rupture volume': [
                    None,
                    2]}},
        'phenotypes_names': [
            'Simple Live',
            'Standard apoptosis model',
            'Standard necrosis model']}}
