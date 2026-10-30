  // the tolerance of the pivot rules of NetworkSimplex: zero for an exact
  // cost type, otherwise relative to the largest of the n costs
  template <typename C>
  inline C ns_pivot_eps(const std::vector<C>& c, int n) {
    if (std::numeric_limits<C>::is_exact) return C(0);
    C m(1);
    for (int i = 0; i != n; ++i) m = std::max(m, c[i] < 0 ? C(-c[i]) : c[i]);
    return m * std::numeric_limits<C>::epsilon() * C(1e6);
  }

  // the tolerance of the flow left on the artificial arcs: zero for an exact
  // flow type, otherwise relative to the largest of the n supplies
  template <typename V>
  inline V ns_flow_eps(const std::vector<V>& s, int n) {
    if (std::numeric_limits<V>::is_exact) return V(0);
    V m(1);
    for (int i = 0; i != n; ++i) m = std::max(m, s[i] < 0 ? V(-s[i]) : s[i]);
    return m * std::numeric_limits<V>::epsilon() * V(1e6);
  }

  // how far the potential of the root of NetworkSimplex may move from zero,
  // which it does when a pivot updates the part of the tree that is not the
  // subtree it moves, given the cost of the artificial arcs: half of it for
  // an exact type, which keeps the sums of the pivot rules within the range
  // of the type, the cost itself otherwise, which costs the potentials no
  // more than one bit of precision
  template <typename C>
  inline C ns_pi_lim(C art) {
    return std::numeric_limits<C>::is_exact ? art / 2 : art;
  }

  // the n + 1 potentials of NetworkSimplex shifted so that the one of the
  // root, the last, is zero
  template <typename C>
  inline void ns_pi_align(std::vector<C>& pi, int n) {
    const C r = pi[n];
    if (r == C(0)) return;
    for (int i = 0; i != n; ++i) pi[i] -= r;
    pi[n] = C(0);
  }
