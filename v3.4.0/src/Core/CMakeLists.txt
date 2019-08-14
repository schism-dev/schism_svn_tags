
### Algorithm options. Everything done using make_define_options will be applied globally.
make_define_option("USE_HA"  "Use the Harmonic Analsis option" OFF "harm.F90")

add_library(core ${local_extra_code}
                 elfe_glbl.F90
                 elfe_msgp.F90
                 selfe_assert.F90
           )
list( APPEND locallibs core)
set(selfemodlibs ${locallibs} ${selfemodlibs} PARENT_SCOPE)




