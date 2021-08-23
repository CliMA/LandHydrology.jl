name: JuliaFormatter

on:
  push:
    branches:
      - main
      - trying
      - staging
    tags: '*'
  pull_request:

jobs:
  format:
    runs-on: ubuntu-latest
    timeout-minutes: 30
    steps:
    - name: Cancel Previous Runs
      uses: styfle/cancel-workflow-action@0.4.0
      with:
        access_token: ${{ github.token }}

    - uses: actions/checkout@v2.2.0

    - uses: dorny/paths-filter@v2.9.1
      id: filter
      with:
        filters: |
          julia_file_change:
            - added|modified: '**.jl'
    - uses: julia-actions/setup-julia@latest
      if: steps.filter.outputs.julia_file_change == 'true'
      with:
        version: 1.6.0

    - name: Apply JuliaFormatter
      if: steps.filter.outputs.julia_file_change == 'true'
      run: |
        julia --project=.dev .dev/climaformat.jl .
    - name: Check formatting diff
      if: steps.filter.outputs.julia_file_change == 'true'
      run: |
        git diff --color=always --exit-code
