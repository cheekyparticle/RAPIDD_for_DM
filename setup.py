from setuptools import setup, find_packages

setup(
    name='rapidd',  # Replace with your package's name
    version='1.0.0',  # Replace with your package's version
    packages=['rapidd'],  # Automatically find all packages in the directory
    description='A code for Direct Detection',  # Replace with your description
    long_description=open('README.md').read(),  # Include a detailed description from the README
    long_description_content_type='text/markdown',  # Markdown is supported in the setup metadata
    author='Andrew Cheek',  # Replace with your name
    author_email='acheek@sjtu.edu.cn',  # Replace with your email
    license='MIT',  # Replace with your chosen license
    url='https://github.com/cheekyparticle/RAPIDD_for_DM',  # Replace with your package link
    package_data={
        'rapidd': ['../lib/build/libRAPIDD.so'],
    },
)
