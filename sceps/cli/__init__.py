"""Command-line entry points for scEPS.

Each module here exposes a ``main()`` that is wired to a console script in
``pyproject.toml``:

==============================  ====================================
command                         module
==============================  ====================================
``sceps``                       :mod:`sceps.cli.run`
``sceps-cluster-neighborhood``  :mod:`sceps.cli.cluster_neighborhood`
``sceps-aggregate``             :mod:`sceps.cli.aggregate`
``sceps-corr``                  :mod:`sceps.cli.corr`
``sceps-generate-test-data``    :mod:`sceps.cli.generate_test_data`
==============================  ====================================
"""
