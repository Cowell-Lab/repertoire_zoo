from connexion import AsyncApp

import repertoire_zoo.bird as br

def post_greeting(name: str):  # Paramaeters are automatically unpacked
    return f"Hello {name}\n", 200          # Responses are automatically serialized

def gene_usage(body):
    print(body)
    return { "image_path": "123.png" }, 200, {"Content-Type": "application/json"}

app = AsyncApp(__name__)
app.add_api("/repertoire_zoo/openapi/vdjserver-image.yaml")
