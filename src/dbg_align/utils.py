def display_mermaid_in_jupyter(mermaid_content: str):
    mm(mermaid_content)

# method to display mermaid in ipynb https://gist.github.com/MLKrisJohnson/2d2df47879ee6afd3be9d6788241fe99
import base64
from IPython.display import Image, display

def mm_ink(graphbytes):
  """Given a bytes object holding a Mermaid-format graph, return a URL that will generate the image."""
  base64_bytes = base64.b64encode(graphbytes)
  base64_string = base64_bytes.decode("ascii")
  return "https://mermaid.ink/img/" + base64_string

def mm_display(graphbytes):
  """Given a bytes object holding a Mermaid-format graph, display it."""
  display(Image(url=mm_ink(graphbytes)))

def mm(graph):
  """Given a string containing a Mermaid-format graph, display it."""
  graphbytes = graph.encode("ascii")
  mm_display(graphbytes)

import graphviz

def display_graphviz(graph):
    """Render and display a Graphviz graph within a Jupyter Notebook."""
    # Render the graph to a file (SVG or PNG can be used here)
    # Note: You might need to adjust the directory path or ensure it exists
    filename = graph.render(filename='temp_graph', format='png', cleanup=True)
    # Display the image in the notebook
    display(Image(filename))