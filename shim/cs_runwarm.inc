
  public:

    /// re-optimizes from the flow and the potentials of the last run
    /** Re-optimizes from the flow and the potentials of the last run instead
     * of starting from a flow found anew and the zero potentials: the flow
     * of the last run is put within the new bounds, and what this and a
     * change of the supplies leave at the nodes is sent where it is due
     * [see feasibleFlow()], which after a change of the costs alone is
     * nothing; the potentials are scaled back to the costs the algorithm
     * works on, and the scaling starts from the largest amount by which a
     * residual arc violates the optimality conditions under them, rather
     * than from the largest cost (the first phase being the one after it,
     * since the pair is already that close to optimal), so that a small
     * change costs a few phases at small epsilon. It does what run() does
     * if the last run did not end at an optimal flow, with lower bounds,
     * with an infinite upper bound, which init() bounds by the deficits,
     * when the supplies do not sum to zero and when no flow meets the new
     * bounds and supplies, which run() then proves. */

    ProblemType runWarm(Method method = PARTIAL_AUGMENT, int factor = 16) {
      if (!_warm || _has_lower || _sum_supply != 0)
        return runCold(method, factor);

      const Value MAX = std::numeric_limits<Value>::max();
      Value sum = 0;
      for (int i = 0; i != _root; ++i) sum += _supply[i];
      if (sum != 0) return runCold(method, factor);
      for (int j = 0; j != _res_arc_num; ++j)
        if (_forward[j] && _source[j] != _root && _upper[j] >= MAX)
          return runCold(method, factor);
      if (!feasibleFlow(true, Value(0), _sup_node_num))
        return runCold(method, factor);

      _alpha = factor;
      const LargeCost scale = static_cast<LargeCost>(_res_node_num) * _alpha;
      for (int i = 0; i != _res_node_num; ++i) {
        _excess[i] = 0;
        _pi[i] *= scale;
      }

      // the costs the algorithm works on, and how far the potentials are
      // from optimal
      _epsilon = 0;
      for (int j = 0; j != _res_arc_num; ++j) {
        const int u = _source[j], v = _target[j];
        if (u == _root || v == _root) continue;
        _cost[j] = static_cast<LargeCost>(_scost[j]) * scale;
      }
      for (int j = 0; j != _res_arc_num; ++j) {
        if (_res_cap[j] <= 0) continue;
        const LargeCost rc = _cost[j] + _pi[_source[j]] - _pi[_target[j]];
        if (-rc > _epsilon) _epsilon = -rc;
      }

      // the pair is already epsilon-optimal, hence the first phase is the
      // one after epsilon, as the loop of start() would take it, rather than
      // a refinement to where the pair already is
      _epsilon = _epsilon >= _alpha ? _epsilon / _alpha
                                    : (_epsilon >= 1 ? 1 : _epsilon);

      _max_rank = _alpha * _res_node_num;
      _buckets.resize(_max_rank);

      start(method);
      _warm = true;
      return OPTIMAL;
    }

  private:

    // run(), after which the flow and the potentials are there to reuse if
    // it has ended at an optimal flow [see runWarm()]
    ProblemType runCold(Method method, int factor) {
      const ProblemType pt = run(method, factor);
      _warm = (pt == OPTIMAL);
      return pt;
    }

    bool _warm = false;

