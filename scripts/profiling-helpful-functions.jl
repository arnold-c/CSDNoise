tinf = @snoop_inference create_optimization_scenarios(specification_vecs)
staleinstances(tinf)
itrigs = inference_triggers(tinf)
mtrigs = accumulate_by_source(Method, itrigs)
modtrigs = SnoopCompile.parcel(mtrigs)
mtrig = mtrigs[2]
summary(mtrig)
itrig = mtrig.itrigs[1]
suggest(itrig)
ascend(itrig)
edit(itrig)

fg = flamegraph(tinf)
using ProfileView
ProfileView.view(fg)

using JET, Cthulhu
report_callees(inference_triggers(tinf))

@code_warntype create_optimization_scenarios(specification_vecs)

@profview create_optimization_scenarios(specification_vecs)
descend_clicked()
