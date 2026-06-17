# Code formatting

1. Code lines should not be longer than 80 columns.
2. Continutation characters will always be placed in the 80th column.
3. Indents are always 4 spaces.
4. Variables should explain what they are. You should be able to understand what the variable is for at a glance.
5. Function and subroutines names should be similarly self-explanatory
6. Functions and subroutines should contain a doc string.
7. When importing modules, use "use modulename, only: var_1, var_2, var_3"
8. Arguments in a function or subroutine should be ordered as: input variables, output variables.
9. Always use intent.
10. Always use doubles.

# Doc strings
Doc strings are contained within the function/subroutine. Each function should have a
doc string that is indented on the same level as the function/subroutine. Doc strings 
should not include change logs. They should include inputs and outputs (subroutines
should specify what the goals of the changed variables are).

```fortran
function blah
!> Human readable doc string (one line minimum)
!> Inputs: 
!>     variable, type
!> Outputs:
!>     variable, type
    code
end function
```

# Structure formatting

1. Subroutines/functions should be contained within a module.
2. Modules should be grouped into concepts (e.g. all raytracing, reflection spectrum)
3. Modules can be grouped together into broader concepts for file structure (e.g. convolution functions).
4. The principal subroutine should be structured such that no actual computation happens within that subroutine. 
5. interface.f90 exists to provide a wrapper (bind C) for foreign function calling. This also requires reltrans to act as a library.
