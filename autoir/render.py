import os
from jinja2 import Environment, FileSystemLoader
from typing import Dict, Any


def render_input_file(params: Dict[str, Any]) -> str:
    """
    Render the LAMMPS input file based on the provided parameters.

    :param params: Dictionary containing job parameters
    :return: Rendered input file content as a string
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
    template_dir = os.path.join(current_dir, "lammps")
    env = Environment(loader=FileSystemLoader(template_dir))

    template = env.get_template("template.in")

    rendered_content = template.render(job=params)

    return rendered_content


def write_input_file(params: Dict[str, Any], output_path: str) -> None:
    """
    Write the rendered LAMMPS input file to the specified output path.

    :param params: Dictionary containing job parameters
    :param output_path: Path where the rendered file should be written
    """
    rendered_content = render_input_file(params)

    with open(output_path, "w") as f:
        f.write(rendered_content)
