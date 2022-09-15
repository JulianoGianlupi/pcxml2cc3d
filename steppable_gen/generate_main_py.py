from pathlib import Path
from os.path import join, isdir

def _main_py_header():
    return "from cc3d import CompuCellSetup\n#------\n\n"


def _main_py_tail():
    return "\nCompuCellSetup.run()\n"


def _steppable_import(step_file, step_name):
    return f"from {step_file} import {step_name}\n"


def _register_steppable(step_name):
    # do not use frequency when generating the import. It will be set on the steppable file itself
    return f"CompuCellSetup.register_steppable(steppable={step_name}())\n\n"


def generate_main_py_string(step_file, step_names):
    full = _main_py_header()

    for name in step_names:
        imp = _steppable_import(step_file, name)
        reg = _register_steppable(name)

        full += imp + reg

    full += _main_py_tail()

    return full

def write_main_py_file(path, filename, main_string):
    path = Path(path)

    if not isdir(path):
        path.mkdir(parents=True)

    with open(join(path, filename), "w+") as f:
        f.write(main_string)

def generate_main_python(path, filename, step_file, step_names):

    if ".py" in step_file:
        step_file = step_file.replace(".py", "")

    full_string = generate_main_py_string(step_file, step_names)

    write_main_py_file(path, filename, full_string)


if __name__ == "__main__":

    steppables = ["steppableA", "steppableB"]
    generate_main_python("D:/test_pc2cc3d/Simulation", "main_test.py", "steppable_test.py", steppables)
