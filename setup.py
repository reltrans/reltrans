import sys
import os
import pathlib
from setuptools import setup

package_name = "pyreltrans"
version = "0.1.0"

setup(
    # TODO: put some actual contact information here
    author="Reltrans Contributors",
    author_email="todo@example.com",
    python_requires=">=3.6",
    install_requires=["xspectrampoline", "numpy"],
    license="MIT",
    name=package_name,
    version=version,
    packages=[package_name],
    include_package_data=True,
    package_data={package_name: ["libpyreltrans.so"]}
)
