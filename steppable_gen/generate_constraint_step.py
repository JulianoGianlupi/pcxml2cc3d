from gen_functions import generate_steppable, steppable_imports

# TODO: have volume, surface, etc, be initialized to their proper cell property

def cell_type_constraint(ctype, this_type_dicts):

    loop = f"\t\tfor cell in self.cell_list_by_type(self.{ctype.upper()}):\n"
    lines = []
    full = loop
    for cell_dict in this_type_dicts:
        for key, value in cell_dict.items():
            if type(value) == str:
                line = f"\t\t\tcell.dict['{key}']='{value}'\n"
            else:
                line = f"\t\t\tcell.dict['{key}']={value}\n"
            full += line
    return full


def get_dicts_for_type(ctype, cell_dicts):
    ds = []
    for d in cell_dicts:
        ds.append(d[ctype])
    return ds


def generate_constraint_loops(cell_types, cell_dicts):
    loops = "\n"
    for ctype in cell_types:
        this_type_dicts = get_dicts_for_type(ctype, cell_dicts)
        loop = cell_type_constraint(ctype, this_type_dicts)
        loops += loop
    return loops


def generate_constraint_steppable(cell_types, cell_dicts, first=True):
    if first:
        imports = steppable_imports()
    else:
        imports = ""

    loops = generate_constraint_loops(cell_types, cell_dicts)

    constraint_step = generate_steppable("Constraints", 1, False, minimal=True, additional_start=loops)

    return constraint_step


if __name__=="__main__":
    dicts = [{"a": {"exemplo":1, "asdasd":3}, "b": {"kkkkkk": "popop"}},
             {"a": {"exeo":1, "uip":3}, "b": {"oooooo": "lololo"}}]
    print(generate_constraint_steppable(["a", "b"], dicts))


