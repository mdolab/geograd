__version__ = "0.9.0"
from .libgeograd import geograd_parallel, triangles, triangles_db
from .libgeograd import geograd as geograd_serial
from .libgeograd_complex import geograd_parallel as geograd_parallel_complex
from .libgeograd_complex import geograd as geograd_serial_complex
from .libgeograd_complex import triangles as triangles_complex
