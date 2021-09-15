from setuptools import setup, find_packages, Extension


setup(
    name='spk',
    install_requires=['invoke'],
    packages=find_packages(),
    license='GPL-2',

    classifiers=[

        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Natural Language :: English',
        'Programming Langauge :: Python :: 3 :: Only',
        'Programming Langauge :: Python :: Implementation :: CPython'

    ],
    ext_modules=[
        Extension(
            'spk',
            ['spkmeansmodule.c', 'spkmeans.c'],
        ),
    ]
)
