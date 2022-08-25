from setuptools import setup
import re

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open("geograd/__init__.py").read(),
)[0]

setup(
    name="geograd",
    version=__version__,
    description="Geograd is a package for efficiency computing triangulated surface constraints in parallel",
    long_description="""Documentation page TBA
      Citation
      --------
      Please cite geograd in any publication for which you find it useful.
      For more background, theory, and figures, see the [geograd journal article](https://arc.aiaa.org/doi/10.2514/1.J058366).
      B. J. Brelje, Anibal, J. L, Yildirim, A., Mader, C. A., and Martins, J. R. R. A., “Flexible Formulation of Spatial Integration Constraints in Aerodynamic Shape Optimization”, in AIAA Journal, 2020. 
      @article{Brelje2020a,
      author = {Benjamin J. Brelje and Joshua Anibal and Anil Yildirim and Charles A. Mader and Joaquim R. R. A. Martins},
      doi = {10.2514/1.J058366},
      journal = {AIAA Journal},
      month = {June},
      number = {6},
      pages = {2571--2580},
      title = {Flexible Formulation of Spatial Integration Constraints in Aerodynamic Shape Optimization},
      volume = {58},
      year = {2020}
    }
      """,
    long_description_content_type="text/markdown",
    keywords="geometric constraints optimization",
    author="",
    author_email="",
    url="https://github.com/mdolab/geograd",
    license="None",
    packages=[
        "geograd",
    ],
    package_data={"geograd": ["*.so"]},
    install_requires=["numpy>=1.16"],
    extras_require={
        "testing": ["numpy>=1.16", "numpy-stl", "openmdao>=2.1", "mpi4py>=3.0"],
        "docs": ["sphinx-mdolab-theme"],
    },
    classifiers=["Operating System :: Linux", "Programming Language :: Python, Fortran"],
)
