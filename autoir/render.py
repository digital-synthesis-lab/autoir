import os
from jinja2 import Environment, FileSystemLoader
from typing import Dict, Any

def render_input_file(params: Dict[str, Any]) -> str:
    """
    Render the LAMMPS input file based on the provided parameters.
    
    :param params: Dictionary containing job parameters
    :return: Rendered input file content as a string
    """
    # Determine which template to use based on the ensemble
    template_name = "npt.in" if params.get("pressure") is not None else "nvt.in"
    
    # Set up Jinja2 environment
    current_dir = os.path.dirname(os.path.abspath(__file__))
    template_dir = os.path.join(current_dir, "lammps")
    env = Environment(loader=FileSystemLoader(template_dir))
    
    # Load the template
    template = env.get_template(template_name)
    
    # Render the template with the provided parameters
    rendered_content = template.render(job=params)
    
    return rendered_content

def write_input_file(params: Dict[str, Any], output_path: str) -> None:
    """
    Write the rendered LAMMPS input file to the specified output path.
    
    :param params: Dictionary containing job parameters
    :param output_path: Path where the rendered file should be written
    """
    rendered_content = render_input_file(params)
    
    with open(output_path, 'w') as f:
        f.write(rendered_content)

# Example usage:
if __name__ == "__main__":
    # Example parameters
    params = {
        "temperature": 300,
        "pressure": 1.0,  # Remove this line for NVT ensemble
        "seed": 12345,
        "equi_steps": 10000,
        "prod_steps": 100000
    }
    
    # Render and print the content
    print(render_input_file(params))
    
    # Optionally, write to a file
    # write_input_file(params, "output.in")
