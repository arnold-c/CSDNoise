#include "manuscript.typ"

#counter(page).update(1)
#counter(figure.where(kind: table)).update(0)
#counter(figure.where(kind: image)).update(0)
#counter(math.equation).update(0)
#counter(footnote).update(0)

#include "supplemental-appendix.typ"
