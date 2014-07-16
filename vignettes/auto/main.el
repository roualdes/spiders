(TeX-add-style-hook "main"
 (lambda ()
    (LaTeX-add-bibliographies
     "refs")
    (LaTeX-add-labels
     "eq:null"
     "eq:lrt"
     "eq:llComp")
    (TeX-run-style-hooks
     "latex2e"
     "art12"
     "article"
     "12pt"
     "RnwOpts")))

