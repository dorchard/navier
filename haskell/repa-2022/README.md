  * `computeP` places in a monad??? But any monad????? Don't get. Paper lies.
    Using `computeS` instead, but less efficient.
  * Think about where to force (compute), what reps combinators should use (e.g.
    `Params.ignoreBoundary`)
  * Similarly, had to use `foldS` because `foldP` requires a monad. Why doesn't
    the upgrade paper describe the impacts of this change, or how to recover
    original behaviour?
  * I think I have to rewrite the *whole* thing to be parametric in some
    monad...
