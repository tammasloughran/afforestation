import setuptools

with open('requirements.txt', 'r') as requirements_file:
    requirements_list = requirements_file.read().splitlines()

setuptools.setup(
        name='example',
        packages=['example'],
        install_requires=requirements_list)

