# Python interface

## Setup

It is possible to call Julia code from Python using the dedicated package [juliacall](https://cjdoris.github.io/PythonCall.jl/stable/juliacall/). This page is a simple tutorial detailing how to set it up and use it for CrystalNets.

In the folder where you intend to put your python script, add a new file called `juliapkg.json` with the following content:

```json
{
    "julia": "1",
    "packages": {
        "CrystalNets": {
            "uuid": "7952bbbe-a946-4118-bea0-081a0932faa9",
            "version": "0.5"
        }
    }
}
```

then install juliacall with `pip install juliacall`. From there, you can call CrystalNets by starting your python script with the following three lines:

```julia
import juliacall
jl = juliacall.newmodule("TheNameOfYourModule") # put whatever name here
jl.seval("using CrystalNets")
```

!!! note
    The very first time these three lines are run in your system, julia will install CrystalNets and its dependencies. This can take a while, but it will only occur once.

Afterwards, you can call any function of CrystalNets.jl (and Julia in general) by prefixing them with `jl.`, for instance:

```python
>>> jl.determine_topology("/path/to/MIL-53.cif")
[ Error: Atom ?, used in a bond, has either zero or multiple placements in the CIF file. This invalidates all bonds from the file, which will thus be discarded.
[ Warning: Guessing bonds with custom algorithm (from Chemfiles and VMD). This may take a while for big structures and may be inexact.
[ Info: To avoid guessing bonds, use a file format that contains the bonds.
Export of input is enabled: saving file at /tmp/input_MIL-53.vtf
Export of subnet_Auto is enabled: saving file at /tmp/subnet_Auto_MIL-53.vtf
Julia: rna
```

Compare this to the equivalent REPL execution in Julia:

```python
julia> using CrystalNets

julia> determine_topology("/path/to/MIL-53.cif")
[ Error: Atom ?, used in a bond, has either zero or multiple placements in the CIF file. This invalidates all bonds from the file, which will thus be discarded.
[ Warning: Guessing bonds with custom algorithm (from Chemfiles and VMD). This may take a while for big structures and may be inexact.
[ Info: To avoid guessing bonds, use a file format that contains the bonds.
Export of input is enabled: saving file at /tmp/input_MIL-53.vtf
Export of subnet_Auto is enabled: saving file at /tmp/subnet_Auto_MIL-53.vtf
rna
```

The same warnings are printed at the beginning, followed by the same exports. The only difference is that the result, **rna** (the topology of MIL-53), is prefixed by `Julia: ` in the Python output.

!!! tip
    You can disable the warnings and/or the exports by adding the following lines after the three initial lines above:
    ```python
    jl.CrystalNets.toggle_warning(False) # to disable warnings
    jl.CrystalNets.toggle_export(False) # to disable exports
    ```

## Usage

Let's now consider a programmatic use-case where the goal is to identify the topology of a complex MOF structure according the [`SingleNodes`](@ref Clustering) and [`AllNodes`](@ref Clustering) clusterings. The main structure may contain interpenetrating substructures.

The function is expected to error if the topologies are different between the two clusterings. Otherwise, it returns a list of pairs whose first element is the dimensionality of the subnet and the second element is the name of the corresponding topology. If there is no known name, the topological genome is used instead.

The Python code is the following:
```python
def identify_topology(cif):
    """Return a list of pairs (dimensionality, topology) for each substructure of the file"""
    options = jl.CrystalNets.Options(structure=jl.StructureType.MOF)
    # Since the structure is specified as a MOF, the default clusterings are AllNodes and SingleNodes
    result = jl.determine_topology(cif, options) # Main call
    # for each x in result:
    # * x[0] is the topology of the substructure.
    # * x[1] is the catenation multiplicity of this subnet.
    return [check_unique_topology(x[0]) for x in result]

def check_unique_topology(result):
    singlenodes = result[jl.Clustering.SingleNodes] # topology for SingleNodes
    allnodes = result[jl.Clustering.AllNodes] # topology for AllNodes
    if singlenodes != allnodes:
        raise Exception("SingleNodes result "+str(singlenodes)+" != AllNodes result "+str(allnodes))
    return (jl.ndims(singlenodes.genome), str(singlenodes))
```

Let's try this on a few examples:

- HKUST-1 has a 3-dimensional topology **tbo**:
  ```python
  >>> identify_topology("/path/to/HKUST-1.cif")
  [(3, 'tbo')]
  ```
- [ODIZIK](https://dx.doi.org/10.5517/cc63w16) is composed of two intepenetrating 2-dimensional topologies **sql**
  ```python
  >>> identify_topology("/path/to/ODIZIK.cif")
  [(2, 'sql'), (2, 'sql')]
  ```
- MIL-53 is has topology **bpq** with `SingleNodes` clustering and **rna** with `AllNodes` so the function errors:
  ```python
  >>> identify_topology("/path/to/MIL-53.cif")
  Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "<stdin>", line 11, in identify_topology
  File "<stdin>", line 5, in check_unique_topology
  Exception: SingleNodes result bpq != AllNodes result rna
  ```
