Step 1 : Go to the directory with the files (fem.jl, manifest.toml .....)

**julia --project=.**


Step 2: Now run these within the REPL. These two lines of code are not required after the first time

**import Pkg**
**Pkg.instantiate()**


Step 3 : Run the desired file in the REPL 

**include("fem.jl")**
**main()**

or

**include("mpm.jl")**
**main()**





