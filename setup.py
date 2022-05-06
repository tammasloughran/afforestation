import setuptools

with open('requirements.txt', 'r') as requirements_file:
    requirements_list = requirements_file.read().splitlines()

setuptools.setup(
        name='analysis',
        packages=['analysis'],
        install_requires=requirements_list)

